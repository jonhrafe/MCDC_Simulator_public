//!  GPU simulation (Phase 1: free diffusion + PGSE signal) ====================
/*! \file gpu_simulation.cu
 *  CUDA implementation behind gpu_simulation.cuh. This translation unit is
 *  Eigen-free and self-contained (only POD in/out). It replicates MCDC's exact
 *  internal-unit phase math (scheme used as-stored, giro in rad/(ms*T), positions
 *  in mm, times in ms) so the GPU result matches the validated CPU free-diffusion
 *  signal. Structure mirrors CUDASynth's deviceDataSynth; the numbers are MCDC's.
 *==============================================================================*/
#include "gpu_simulation.cuh"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath>
#include <cstdio>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//! Gyromagnetic ratio, rad/(ms*T) -- identical to MCDC/CUDASynth constants.h.
__device__ __constant__ float GIRO = 267.51525e3f;

#define CUDA_CHECK(call)                                                         \
    do {                                                                         \
        cudaError_t _e = (call);                                                 \
        if (_e != cudaSuccess) {                                                 \
            std::fprintf(stderr, "[GPU][ERROR] %s:%d %s\n", __FILE__, __LINE__,  \
                         cudaGetErrorString(_e));                                \
            return int(_e);                                                      \
        }                                                                        \
    } while (0)

//! PGSE gradient impulse for a timestep [tLast,t]. Verbatim port of MCDC's
//! PGSESequence::getGradImpulse (pgsesequence.cpp) -- scheme values used directly,
//! no unit conversion (MCDC stores/uses internal units).
__device__ void getGradImpulse(float t, float tLast, const float g[3], float G,
                               float Delta, float delta, float te, float Gdt[3])
{
    Gdt[0] = Gdt[1] = Gdt[2] = 0.0f;

    if (!(t >= 0.0f && t <= te)) return;

    float pad = (te - Delta - delta) / 2.0f;
    if ((t < pad) || (t > te - pad)) return;

    float firstBlockStart = pad;
    float firstBlockEnd   = pad + delta;

    float sgn = 1.0f;
    if (t >= firstBlockStart && t < firstBlockEnd) {
        if (tLast < firstBlockStart) {
            float dt = t - firstBlockStart;
            for (int j = 0; j < 3; ++j) Gdt[j] = sgn * G * dt * g[j];
            return;
        }
    }

    float secondBlockStart = pad + Delta;
    float secondBlockEnd   = pad + Delta + delta;

    if (t >= firstBlockEnd && t < secondBlockStart) {
        if (tLast < firstBlockEnd) {
            float dt = firstBlockEnd - tLast;
            for (int j = 0; j < 3; ++j) Gdt[j] = sgn * G * dt * g[j];
            return;
        }
        return;
    }

    if (t >= secondBlockStart) sgn = -1.0f;

    if (t >= secondBlockEnd) {
        if (tLast < secondBlockEnd) {
            float dt = secondBlockEnd - tLast;
            for (int j = 0; j < 3; ++j) Gdt[j] = sgn * G * dt * g[j];
            return;
        }
        return;
    }

    if ((t >= secondBlockStart) && (tLast < secondBlockStart)) {
        float dt = t - secondBlockStart;
        for (int j = 0; j < 3; ++j) Gdt[j] = sgn * G * dt * g[j];
        return;
    }

    for (int j = 0; j < 3; ++j) Gdt[j] = sgn * G * (t - tLast) * g[j];
}

//! One thread per walker: free random walk (fixed step length, random direction),
//! writing the trajectory [walker][step][xyz] = (T+1)*3 floats.
__global__ void freeDiffusionKernel(float* traj, unsigned N, unsigned T,
                                    float step_length, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);              // counter-based, per-walker subsequence

    float px = 0.0f, py = 0.0f, pz = 0.0f;
    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = tr[1] = tr[2] = 0.0f;

    for (unsigned t = 1; t <= T; ++t) {
        // isotropic direction: 3 normals normalized to the unit sphere
        float dx = curand_normal(&st);
        float dy = curand_normal(&st);
        float dz = curand_normal(&st);
        float inv = rsqrtf(dx*dx + dy*dy + dz*dz + 1e-30f);
        px += step_length * dx * inv;
        py += step_length * dy * inv;
        pz += step_length * dz * inv;
        float* o = tr + (size_t)t * 3;
        o[0] = px; o[1] = py; o[2] = pz;
    }
}

//! Elastic reflection of a unit direction about a unit surface normal: d -= 2(d.n)n.
__device__ inline void reflectDir(float d[3], const float n[3])
{
    float dn = d[0]*n[0] + d[1]*n[1] + d[2]*n[2];
    d[0] -= 2.0f*dn*n[0];
    d[1] -= 2.0f*dn*n[1];
    d[2] -= 2.0f*dn*n[2];
}

//! Phase 2a/2c: one thread per walker, restricted diffusion among impermeable
//! spheres (spheres[i*4+0..3] = cx,cy,cz,r). Walker starts uniformly inside a
//! random sphere; each step takes the NEAREST surface crossing over all spheres
//! and reflects (multi-bounce). Brute-force over all spheres (no grid yet).
__global__ void spheresKernel(float* traj, unsigned N, unsigned T, float step_length,
                              const float* spheres, unsigned n_sph, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    // Start uniformly inside a randomly chosen sphere (equal-radius packings ->
    // uniform over the intra volume).
    unsigned si = curand(&st) % n_sph;
    float cx = spheres[si*4+0], cy = spheres[si*4+1], cz = spheres[si*4+2], r0 = spheres[si*4+3];
    float ox, oy, oz;
    do {
        ox = (2.0f*curand_uniform(&st) - 1.0f) * r0;
        oy = (2.0f*curand_uniform(&st) - 1.0f) * r0;
        oz = (2.0f*curand_uniform(&st) - 1.0f) * r0;
    } while (ox*ox + oy*oy + oz*oz > r0*r0);
    float px = cx + ox, py = cy + oy, pz = cz + oz;

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = px; tr[1] = py; tr[2] = pz;

    const float EPS = 1e-6f * r0;
    for (unsigned t = 1; t <= T; ++t) {
        float dx = curand_normal(&st), dy = curand_normal(&st), dz = curand_normal(&st);
        float inv = rsqrtf(dx*dx + dy*dy + dz*dz + 1e-30f);
        float dir[3] = { dx*inv, dy*inv, dz*inv };
        float remaining = step_length;

        for (int b = 0; b < 16 && remaining > 0.0f; ++b) {
            float best_t = remaining; int best = -1;
            for (unsigned i = 0; i < n_sph; ++i) {
                float ux = px - spheres[i*4+0], uy = py - spheres[i*4+1], uz = pz - spheres[i*4+2];
                float ri = spheres[i*4+3];
                float bq = ux*dir[0] + uy*dir[1] + uz*dir[2];
                float cq = ux*ux + uy*uy + uz*uz - ri*ri;
                float discr = bq*bq - cq;
                if (discr > 0.0f) {
                    float sq = sqrtf(discr);
                    float t1 = -bq - sq, t2 = -bq + sq;       // first crossing > EPS
                    float tc = (t1 > EPS) ? t1 : ((t2 > EPS) ? t2 : -1.0f);
                    if (tc > EPS && tc < best_t) { best_t = tc; best = int(i); }
                }
            }
            if (best < 0) {
                px += remaining*dir[0]; py += remaining*dir[1]; pz += remaining*dir[2];
                remaining = 0.0f;
            } else {
                px += best_t*dir[0]; py += best_t*dir[1]; pz += best_t*dir[2];
                float nx = px - spheres[best*4+0], ny = py - spheres[best*4+1], nz = pz - spheres[best*4+2];
                float nn = rsqrtf(nx*nx + ny*ny + nz*nz + 1e-30f);
                float n[3] = { nx*nn, ny*nn, nz*nn };
                reflectDir(dir, n);
                // nudge along the reflected direction (into the valid region;
                // works whether we left our own sphere or bounced off a neighbour)
                px += EPS*dir[0]; py += EPS*dir[1]; pz += EPS*dir[2];
                remaining -= best_t;
            }
        }
        float* o = tr + (size_t)t * 3;
        o[0] = px; o[1] = py; o[2] = pz;
    }
}

//! Phase 2b/2c: one thread per walker, restricted diffusion among impermeable
//! infinite cylinders (cyls[i*7+0..6] = px,py,pz, dx,dy,dz(unit), r). Walker starts
//! in a random cylinder's disk; each step takes the NEAREST radial crossing over
//! all cylinders and reflects off that radial normal (axial motion is free).
__global__ void cylindersKernel(float* traj, unsigned N, unsigned T, float step_length,
                                const float* cyls, unsigned n_cyl, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    // Start in a random cylinder's perpendicular disk.
    unsigned ci = curand(&st) % n_cyl;
    float p0x = cyls[ci*7+0], p0y = cyls[ci*7+1], p0z = cyls[ci*7+2];
    float dx = cyls[ci*7+3], dy = cyls[ci*7+4], dz = cyls[ci*7+5], r0 = cyls[ci*7+6];
    float ax = (fabsf(dx) < 0.9f) ? 1.0f : 0.0f;
    float ay = (fabsf(dx) < 0.9f) ? 0.0f : 1.0f;
    float adb = ax*dx + ay*dy;
    float e1[3] = { ax - adb*dx, ay - adb*dy, -adb*dz };
    float i1 = rsqrtf(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2] + 1e-30f);
    e1[0]*=i1; e1[1]*=i1; e1[2]*=i1;
    float e2[3] = { dy*e1[2]-dz*e1[1], dz*e1[0]-dx*e1[2], dx*e1[1]-dy*e1[0] };
    float rr = r0 * sqrtf(curand_uniform(&st));
    float th = 6.2831853f * curand_uniform(&st);
    float cs = cosf(th)*rr, sn = sinf(th)*rr;
    float px = p0x + cs*e1[0] + sn*e2[0];
    float py = p0y + cs*e1[1] + sn*e2[1];
    float pz = p0z + cs*e1[2] + sn*e2[2];

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = px; tr[1] = py; tr[2] = pz;

    const float EPS = 1e-6f * r0;
    for (unsigned t = 1; t <= T; ++t) {
        float gx = curand_normal(&st), gy = curand_normal(&st), gz = curand_normal(&st);
        float inv = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-30f);
        float dir[3] = { gx*inv, gy*inv, gz*inv };
        float remaining = step_length;

        for (int b = 0; b < 16 && remaining > 0.0f; ++b) {
            float best_t = remaining; int best = -1;
            for (unsigned i = 0; i < n_cyl; ++i) {
                float ex = px - cyls[i*7+0], ey = py - cyls[i*7+1], ez = pz - cyls[i*7+2];
                float Dx = cyls[i*7+3], Dy = cyls[i*7+4], Dz = cyls[i*7+5], ri = cyls[i*7+6];
                float md = ex*Dx + ey*Dy + ez*Dz;
                float nd = dir[0]*Dx + dir[1]*Dy + dir[2]*Dz;
                float a  = 1.0f - nd*nd;
                if (a <= 1e-8f) continue;                 // parallel to this axis
                float cq = (ex*ex + ey*ey + ez*ez - ri*ri) - md*md;
                float bq = (ex*dir[0] + ey*dir[1] + ez*dir[2]) - nd*md;
                float discr = bq*bq - a*cq;
                if (discr > 0.0f) {
                    float sq = sqrtf(discr);
                    float t1 = (-bq - sq) / a, t2 = (-bq + sq) / a;
                    float tc = (t1 > EPS) ? t1 : ((t2 > EPS) ? t2 : -1.0f);
                    if (tc > EPS && tc < best_t) { best_t = tc; best = int(i); }
                }
            }
            if (best < 0) {
                px += remaining*dir[0]; py += remaining*dir[1]; pz += remaining*dir[2];
                remaining = 0.0f;
            } else {
                px += best_t*dir[0]; py += best_t*dir[1]; pz += best_t*dir[2];
                float Dx = cyls[best*7+3], Dy = cyls[best*7+4], Dz = cyls[best*7+5];
                float rx = px - cyls[best*7+0], ry = py - cyls[best*7+1], rz = pz - cyls[best*7+2];
                float rad = rx*Dx + ry*Dy + rz*Dz;
                float n[3] = { rx - rad*Dx, ry - rad*Dy, rz - rad*Dz };
                float nn = rsqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] + 1e-30f);
                n[0]*=nn; n[1]*=nn; n[2]*=nn;
                reflectDir(dir, n);
                px += EPS*dir[0]; py += EPS*dir[1]; pz += EPS*dir[2];
                remaining -= best_t;
            }
        }
        float* o = tr + (size_t)t * 3;
        o[0] = px; o[1] = py; o[2] = pz;
    }
}

//! One block per walker, num_rep measurements strided across threads. Reads the
//! on-device trajectory and accumulates the (T2=inf) signal sum_w cos(phase).
//! Mirrors CUDASynth::deviceDataSynth but uses the scheme values directly (MCDC).
__global__ void signalKernel(const float* traj, unsigned T, float duration,
                             const float* scheme, unsigned num_rep, float* out_DWI)
{
    unsigned w = blockIdx.x;
    float time_step = duration / float(T);
    float dos_pi = 2.0f * float(M_PI);

    const float* tr = traj + (size_t)w * (T + 1) * 3;
    float x_ini[3] = { tr[0], tr[1], tr[2] };

    for (unsigned s = threadIdx.x; s < num_rep; s += blockDim.x) {
        unsigned k = s * 7;
        float g[3] = { scheme[k + 0], scheme[k + 1], scheme[k + 2] };
        float G     = scheme[k + 3];   // MCDC stores/uses internal units directly
        float Delta = scheme[k + 4];
        float delta = scheme[k + 5];
        float te    = scheme[k + 6];

        float Gdt[3], xt[3];
        float phase = 0.0f;
        for (unsigned t = 1; t <= T; ++t) {
            float dt_last = time_step * float(t - 1);
            float dt      = time_step * float(t);
            const float* xp = tr + (size_t)t * 3;
            xt[0] = xp[0] - x_ini[0];
            xt[1] = xp[1] - x_ini[1];
            xt[2] = xp[2] - x_ini[2];
            getGradImpulse(dt, dt_last, g, G, Delta, delta, te, Gdt);
            float val = GIRO * (Gdt[0]*xt[0] + Gdt[1]*xt[1] + Gdt[2]*xt[2]);
            val   = fmodf(val, dos_pi);
            phase = fmodf(phase + val, dos_pi);
        }
        atomicAdd(&out_DWI[s], cosf(phase));
    }
}

int runGpuFreeDiffusion(const GpuSimParams& p, const float* scheme, double* out_DWI)
{
    const size_t traj_floats = (size_t)p.N * (p.T + 1) * 3;
    const size_t traj_bytes  = traj_floats * sizeof(float);
    const size_t scheme_bytes = (size_t)p.num_rep * 7 * sizeof(float);
    const size_t dwi_bytes    = (size_t)p.num_rep * sizeof(float);

    float *d_traj = nullptr, *d_scheme = nullptr, *d_dwi = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj,   traj_bytes));
    CUDA_CHECK(cudaMalloc(&d_scheme, scheme_bytes));
    CUDA_CHECK(cudaMalloc(&d_dwi,    dwi_bytes));
    CUDA_CHECK(cudaMemcpy(d_scheme, scheme, scheme_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_dwi, 0, dwi_bytes));

    const int bs = 256;
    const int gs = int((p.N + bs - 1) / bs);
    freeDiffusionKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, p.seed);
    CUDA_CHECK(cudaGetLastError());

    const int threads = p.num_rep < 1024u ? int(p.num_rep) : 1024;
    signalKernel<<<p.N, threads>>>(d_traj, p.T, p.duration, d_scheme, p.num_rep, d_dwi);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<float> h_dwi(p.num_rep);
    CUDA_CHECK(cudaMemcpy(h_dwi.data(), d_dwi, dwi_bytes, cudaMemcpyDeviceToHost));
    for (unsigned i = 0; i < p.num_rep; ++i) out_DWI[i] = double(h_dwi[i]);

    cudaFree(d_traj);
    cudaFree(d_scheme);
    cudaFree(d_dwi);
    return 0;
}

// Shared tail: upload scheme, run the signal kernel on the (already-filled) device
// trajectory, copy back the DWI, free. Returns 0 / CUDA error code.
static int gpuSynthAndFree(const GpuSimParams& p, float* d_traj, const float* scheme, double* out_DWI)
{
    const size_t scheme_bytes = (size_t)p.num_rep * 7 * sizeof(float);
    const size_t dwi_bytes    = (size_t)p.num_rep * sizeof(float);
    float *d_scheme = nullptr, *d_dwi = nullptr;
    CUDA_CHECK(cudaMalloc(&d_scheme, scheme_bytes));
    CUDA_CHECK(cudaMalloc(&d_dwi,    dwi_bytes));
    CUDA_CHECK(cudaMemcpy(d_scheme, scheme, scheme_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_dwi, 0, dwi_bytes));

    const int threads = p.num_rep < 1024u ? int(p.num_rep) : 1024;
    signalKernel<<<p.N, threads>>>(d_traj, p.T, p.duration, d_scheme, p.num_rep, d_dwi);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<float> h_dwi(p.num_rep);
    CUDA_CHECK(cudaMemcpy(h_dwi.data(), d_dwi, dwi_bytes, cudaMemcpyDeviceToHost));
    for (unsigned i = 0; i < p.num_rep; ++i) out_DWI[i] = double(h_dwi[i]);

    cudaFree(d_scheme);
    cudaFree(d_dwi);
    return 0;
}

int runGpuSpheres(const GpuSimParams& p, const float* spheres, unsigned n_sph,
                  const float* scheme, double* out_DWI)
{
    float *d_traj = nullptr, *d_sph = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj, (size_t)p.N * (p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_sph,  (size_t)n_sph * 4 * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_sph, spheres, (size_t)n_sph * 4 * sizeof(float), cudaMemcpyHostToDevice));

    const int bs = 256, gs = int((p.N + bs - 1) / bs);
    spheresKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, d_sph, n_sph, p.seed);
    CUDA_CHECK(cudaGetLastError());

    int rc = gpuSynthAndFree(p, d_traj, scheme, out_DWI);
    cudaFree(d_traj);
    cudaFree(d_sph);
    return rc;
}

int runGpuCylinders(const GpuSimParams& p, const float* cyls, unsigned n_cyl,
                    const float* scheme, double* out_DWI)
{
    float *d_traj = nullptr, *d_cyl = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj, (size_t)p.N * (p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_cyl,  (size_t)n_cyl * 7 * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_cyl, cyls, (size_t)n_cyl * 7 * sizeof(float), cudaMemcpyHostToDevice));

    const int bs = 256, gs = int((p.N + bs - 1) / bs);
    cylindersKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, d_cyl, n_cyl, p.seed);
    CUDA_CHECK(cudaGetLastError());

    int rc = gpuSynthAndFree(p, d_traj, scheme, out_DWI);
    cudaFree(d_traj);
    cudaFree(d_cyl);
    return rc;
}
