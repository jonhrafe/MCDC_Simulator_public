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

//! Phase 2a: one thread per walker, restricted diffusion inside ONE impermeable
//! sphere. Walkers start uniformly inside; each step elastically reflects off the
//! inner surface (multi-bounce), brute-force (single obstacle, no grid).
__global__ void sphereKernel(float* traj, unsigned N, unsigned T, float step_length,
                             float cx, float cy, float cz, float r, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    // Uniform start inside the sphere (rejection in the bounding cube).
    float px, py, pz;
    do {
        px = (2.0f*curand_uniform(&st) - 1.0f) * r;
        py = (2.0f*curand_uniform(&st) - 1.0f) * r;
        pz = (2.0f*curand_uniform(&st) - 1.0f) * r;
    } while (px*px + py*py + pz*pz > r*r);
    px += cx; py += cy; pz += cz;

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = px; tr[1] = py; tr[2] = pz;

    const float EPS = 1e-6f * r;   // nudge to stay inside after a bounce
    for (unsigned t = 1; t <= T; ++t) {
        float dx = curand_normal(&st), dy = curand_normal(&st), dz = curand_normal(&st);
        float inv = rsqrtf(dx*dx + dy*dy + dz*dz + 1e-30f);
        float dir[3] = { dx*inv, dy*inv, dz*inv };
        float remaining = step_length;

        for (int b = 0; b < 16 && remaining > 0.0f; ++b) {
            float mx = px - cx, my = py - cy, mz = pz - cz;
            float bq = mx*dir[0] + my*dir[1] + mz*dir[2];
            float cq = mx*mx + my*my + mz*mz - r*r;
            float discr = bq*bq - cq;
            float thit = -1.0f;
            if (discr > 0.0f) {
                float t2 = -bq + sqrtf(discr);     // exit root (inside => cq<0 => t2>0)
                if (t2 > EPS && t2 < remaining) thit = t2;
            }
            if (thit < 0.0f) {                     // no surface hit within the step
                px += remaining*dir[0]; py += remaining*dir[1]; pz += remaining*dir[2];
                remaining = 0.0f;
            } else {                               // move to surface, reflect, continue
                px += thit*dir[0]; py += thit*dir[1]; pz += thit*dir[2];
                float nn = rsqrtf((px-cx)*(px-cx) + (py-cy)*(py-cy) + (pz-cz)*(pz-cz) + 1e-30f);
                float n[3] = { (px-cx)*nn, (py-cy)*nn, (pz-cz)*nn };
                reflectDir(dir, n);
                px -= EPS*n[0]; py -= EPS*n[1]; pz -= EPS*n[2];   // nudge inside
                remaining -= thit;
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

int runGpuSphere(const GpuSimParams& p, const GpuSphere& s, const float* scheme, double* out_DWI)
{
    const size_t traj_bytes   = (size_t)p.N * (p.T + 1) * 3 * sizeof(float);
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
    sphereKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, s.cx, s.cy, s.cz, s.r, p.seed);
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
