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
// Free diffusion inside the voxel with boundary conditions. Two positions:
//   pr = real/unwrapped position (the trajectory used for the signal)
//   pv = voxel/wrapped position  (kept inside [vmin,vmax], used for the geometry)
// periodic: pv wraps to the opposite face, pr continues unwrapped (-> free signal,
//           and the mechanism that tiles obstacles in later phases).
// mirror  : the step reflects off the face, so pr is bounded to the box. NOTE this
//           is a folded (physically bounded) convention; the CPU uses an unfolded
//           boundary_mirror scheme, so GPU mirror does NOT yet match CPU mirror.
//           Only periodic BC is validated against the CPU so far.
__global__ void freeDiffusionKernel(float* traj, unsigned N, unsigned T, float step_length,
                                    float vminx, float vminy, float vminz,
                                    float vmaxx, float vmaxy, float vmaxz,
                                    int periodic, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    float vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    // start uniformly in the voxel; pr and pv coincide at t=0
    float pv[3] = { vmin[0] + curand_uniform(&st)*(vmax[0]-vmin[0]),
                    vmin[1] + curand_uniform(&st)*(vmax[1]-vmin[1]),
                    vmin[2] + curand_uniform(&st)*(vmax[2]-vmin[2]) };
    float pr[3] = { pv[0], pv[1], pv[2] };

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];

    const float EPS = 1e-5f * step_length;
    for (unsigned t = 1; t <= T; ++t) {
        float gx = curand_normal(&st), gy = curand_normal(&st), gz = curand_normal(&st);
        float inv = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-30f);
        float dir[3] = { gx*inv, gy*inv, gz*inv };
        float remaining = step_length;

        for (int b = 0; b < 16 && remaining > 0.0f; ++b) {
            // nearest voxel face along dir, from the wrapped position pv
            float best_t = remaining; int axis = -1;
            for (int i = 0; i < 3; ++i) {
                float tw = -1.0f;
                if (dir[i] >  1e-20f) tw = (vmax[i] - pv[i]) / dir[i];
                else if (dir[i] < -1e-20f) tw = (vmin[i] - pv[i]) / dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; axis = i; }
            }
            if (axis < 0) {
                for (int i = 0; i < 3; ++i) { pr[i] += remaining*dir[i]; pv[i] += remaining*dir[i]; }
                remaining = 0.0f;
            } else {
                for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]; pv[i] += best_t*dir[i]; }
                remaining -= best_t;
                if (periodic) {
                    pv[axis] = (dir[axis] > 0.0f) ? vmin[axis] : vmax[axis];   // wrap pv only
                } else {
                    dir[axis] = -dir[axis];                                    // mirror: reflect pr & pv
                    pr[axis] += EPS*dir[axis]; pv[axis] += EPS*dir[axis];      // nudge off the wall
                }
            }
        }
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
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
__device__ inline bool insideAnySphere(float x, float y, float z, const float* spheres, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        float dx = x - spheres[i*4+0], dy = y - spheres[i*4+1], dz = z - spheres[i*4+2];
        float r = spheres[i*4+3];
        if (dx*dx + dy*dy + dz*dz < r*r) return true;
    }
    return false;
}

//! Phase 2a/2c/2e: restricted diffusion among impermeable spheres inside a voxel.
//! Two positions: pv (wrapped, for collision) and pr (unwrapped, for the signal).
//! Each step takes the nearest of {sphere surfaces, voxel faces}; spheres reflect,
//! voxel faces wrap pv (periodic) or reflect (mirror, folded). placement: 0=uniform
//! in voxel, 1=intra (random sphere), 2=extra (voxel, outside all spheres).
__global__ void spheresKernel(float* traj, unsigned N, unsigned T, float step_length,
                              const float* spheres, unsigned n_sph,
                              float vminx, float vminy, float vminz,
                              float vmaxx, float vmaxy, float vmaxz,
                              int periodic, int placement, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    float vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    float pv[3];
    if (placement == 1) {                          // intra: uniform in a random sphere
        unsigned si = curand(&st) % n_sph;
        float cx = spheres[si*4+0], cy = spheres[si*4+1], cz = spheres[si*4+2], r0 = spheres[si*4+3];
        float ox, oy, oz;
        do { ox = (2.0f*curand_uniform(&st)-1.0f)*r0; oy = (2.0f*curand_uniform(&st)-1.0f)*r0; oz = (2.0f*curand_uniform(&st)-1.0f)*r0; }
        while (ox*ox + oy*oy + oz*oz > r0*r0);
        pv[0] = cx + ox; pv[1] = cy + oy; pv[2] = cz + oz;
    } else {                                       // uniform (0) or extra (2): in the voxel
        for (int tries = 0; tries < 10000; ++tries) {
            pv[0] = vmin[0] + curand_uniform(&st)*(vmax[0]-vmin[0]);
            pv[1] = vmin[1] + curand_uniform(&st)*(vmax[1]-vmin[1]);
            pv[2] = vmin[2] + curand_uniform(&st)*(vmax[2]-vmin[2]);
            if (placement != 2 || !insideAnySphere(pv[0], pv[1], pv[2], spheres, n_sph)) break;
        }
    }
    float pr[3] = { pv[0], pv[1], pv[2] };

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];

    const float EPS = 1e-5f * step_length;
    for (unsigned t = 1; t <= T; ++t) {
        float dx = curand_normal(&st), dy = curand_normal(&st), dz = curand_normal(&st);
        float inv = rsqrtf(dx*dx + dy*dy + dz*dz + 1e-30f);
        float dir[3] = { dx*inv, dy*inv, dz*inv };
        float remaining = step_length;

        for (int b = 0; b < 24 && remaining > 0.0f; ++b) {
            float best_t = remaining; int best_sph = -1, best_axis = -1;
            for (unsigned i = 0; i < n_sph; ++i) {     // sphere surfaces (collide in pv)
                float ux = pv[0]-spheres[i*4+0], uy = pv[1]-spheres[i*4+1], uz = pv[2]-spheres[i*4+2];
                float ri = spheres[i*4+3];
                float bq = ux*dir[0] + uy*dir[1] + uz*dir[2];
                float cq = ux*ux + uy*uy + uz*uz - ri*ri;
                float discr = bq*bq - cq;
                if (discr > 0.0f) {
                    float sq = sqrtf(discr);
                    float t1 = -bq - sq, t2 = -bq + sq;
                    float tc = (t1 > EPS) ? t1 : ((t2 > EPS) ? t2 : -1.0f);
                    if (tc > EPS && tc < best_t) { best_t = tc; best_sph = int(i); best_axis = -1; }
                }
            }
            for (int i = 0; i < 3; ++i) {               // voxel faces (in pv)
                float tw = -1.0f;
                if (dir[i] >  1e-20f) tw = (vmax[i]-pv[i])/dir[i];
                else if (dir[i] < -1e-20f) tw = (vmin[i]-pv[i])/dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; best_axis = i; best_sph = -1; }
            }
            // advance both positions to the event
            for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]; pv[i] += best_t*dir[i]; }
            remaining -= best_t;
            if (best_sph >= 0) {                         // reflect off the sphere
                float nx = pv[0]-spheres[best_sph*4+0], ny = pv[1]-spheres[best_sph*4+1], nz = pv[2]-spheres[best_sph*4+2];
                float nn = rsqrtf(nx*nx + ny*ny + nz*nz + 1e-30f);
                float n[3] = { nx*nn, ny*nn, nz*nn };
                reflectDir(dir, n);
                for (int i = 0; i < 3; ++i) { pr[i] += EPS*dir[i]; pv[i] += EPS*dir[i]; }
            } else if (best_axis >= 0) {                 // voxel face
                if (periodic) {
                    pv[best_axis] = (dir[best_axis] > 0.0f) ? vmin[best_axis] : vmax[best_axis];
                } else {
                    dir[best_axis] = -dir[best_axis];
                    pr[best_axis] += EPS*dir[best_axis]; pv[best_axis] += EPS*dir[best_axis];
                }
            } else {
                remaining = 0.0f;                        // no event: full step already taken above
            }
        }
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
    }
}

//! Phase 2b/2c: one thread per walker, restricted diffusion among impermeable
//! infinite cylinders (cyls[i*7+0..6] = px,py,pz, dx,dy,dz(unit), r). Walker starts
//! in a random cylinder's disk; each step takes the NEAREST radial crossing over
//! all cylinders and reflects off that radial normal (axial motion is free).
__device__ inline bool insideAnyCylinder(float x, float y, float z, const float* cyls, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        float ex = x-cyls[i*7+0], ey = y-cyls[i*7+1], ez = z-cyls[i*7+2];
        float Dx = cyls[i*7+3], Dy = cyls[i*7+4], Dz = cyls[i*7+5], r = cyls[i*7+6];
        float md = ex*Dx + ey*Dy + ez*Dz;
        float rad2 = ex*ex + ey*ey + ez*ez - md*md;
        if (rad2 < r*r) return true;
    }
    return false;
}

//! Phase 2b/2c/2e: restricted diffusion among impermeable infinite cylinders inside
//! a voxel. Two positions (pv wrapped / pr unwrapped); each step takes the nearest of
//! {cylinder surfaces, voxel faces}. placement: 0=uniform, 1=intra (random cylinder
//! disk), 2=extra (in the voxel, outside all cylinders).
__global__ void cylindersKernel(float* traj, unsigned N, unsigned T, float step_length,
                                const float* cyls, unsigned n_cyl,
                                float vminx, float vminy, float vminz,
                                float vmaxx, float vmaxy, float vmaxz,
                                int periodic, int placement, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    float vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    float pv[3];
    if (placement == 1) {                          // intra: random cylinder's perpendicular disk
        unsigned ci = curand(&st) % n_cyl;
        float p0x = cyls[ci*7+0], p0y = cyls[ci*7+1], p0z = cyls[ci*7+2];
        float dx = cyls[ci*7+3], dy = cyls[ci*7+4], dz = cyls[ci*7+5], r0 = cyls[ci*7+6];
        float ax = (fabsf(dx) < 0.9f) ? 1.0f : 0.0f, ay = (fabsf(dx) < 0.9f) ? 0.0f : 1.0f;
        float adb = ax*dx + ay*dy;
        float e1[3] = { ax - adb*dx, ay - adb*dy, -adb*dz };
        float i1 = rsqrtf(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2] + 1e-30f);
        e1[0]*=i1; e1[1]*=i1; e1[2]*=i1;
        float e2[3] = { dy*e1[2]-dz*e1[1], dz*e1[0]-dx*e1[2], dx*e1[1]-dy*e1[0] };
        float rr = r0 * sqrtf(curand_uniform(&st)), th = 6.2831853f * curand_uniform(&st);
        float cs = cosf(th)*rr, sn = sinf(th)*rr;
        pv[0] = p0x + cs*e1[0] + sn*e2[0]; pv[1] = p0y + cs*e1[1] + sn*e2[1]; pv[2] = p0z + cs*e1[2] + sn*e2[2];
    } else {                                       // uniform (0) or extra (2)
        for (int tries = 0; tries < 10000; ++tries) {
            pv[0] = vmin[0] + curand_uniform(&st)*(vmax[0]-vmin[0]);
            pv[1] = vmin[1] + curand_uniform(&st)*(vmax[1]-vmin[1]);
            pv[2] = vmin[2] + curand_uniform(&st)*(vmax[2]-vmin[2]);
            if (placement != 2 || !insideAnyCylinder(pv[0], pv[1], pv[2], cyls, n_cyl)) break;
        }
    }
    float pr[3] = { pv[0], pv[1], pv[2] };

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];

    const float EPS = 1e-5f * step_length;
    for (unsigned t = 1; t <= T; ++t) {
        float gx = curand_normal(&st), gy = curand_normal(&st), gz = curand_normal(&st);
        float inv = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-30f);
        float dir[3] = { gx*inv, gy*inv, gz*inv };
        float remaining = step_length;

        for (int b = 0; b < 24 && remaining > 0.0f; ++b) {
            float best_t = remaining; int best_cyl = -1, best_axis = -1;
            for (unsigned i = 0; i < n_cyl; ++i) {     // cylinder surfaces (in pv)
                float ex = pv[0]-cyls[i*7+0], ey = pv[1]-cyls[i*7+1], ez = pv[2]-cyls[i*7+2];
                float Dx = cyls[i*7+3], Dy = cyls[i*7+4], Dz = cyls[i*7+5], ri = cyls[i*7+6];
                float md = ex*Dx + ey*Dy + ez*Dz;
                float nd = dir[0]*Dx + dir[1]*Dy + dir[2]*Dz;
                float a  = 1.0f - nd*nd;
                if (a <= 1e-8f) continue;
                float cq = (ex*ex + ey*ey + ez*ez - ri*ri) - md*md;
                float bq = (ex*dir[0] + ey*dir[1] + ez*dir[2]) - nd*md;
                float discr = bq*bq - a*cq;
                if (discr > 0.0f) {
                    float sq = sqrtf(discr);
                    float t1 = (-bq - sq) / a, t2 = (-bq + sq) / a;
                    float tc = (t1 > EPS) ? t1 : ((t2 > EPS) ? t2 : -1.0f);
                    if (tc > EPS && tc < best_t) { best_t = tc; best_cyl = int(i); best_axis = -1; }
                }
            }
            for (int i = 0; i < 3; ++i) {               // voxel faces (in pv)
                float tw = -1.0f;
                if (dir[i] >  1e-20f) tw = (vmax[i]-pv[i])/dir[i];
                else if (dir[i] < -1e-20f) tw = (vmin[i]-pv[i])/dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; best_axis = i; best_cyl = -1; }
            }
            for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]; pv[i] += best_t*dir[i]; }
            remaining -= best_t;
            if (best_cyl >= 0) {                         // reflect off the cylinder
                float Dx = cyls[best_cyl*7+3], Dy = cyls[best_cyl*7+4], Dz = cyls[best_cyl*7+5];
                float rx = pv[0]-cyls[best_cyl*7+0], ry = pv[1]-cyls[best_cyl*7+1], rz = pv[2]-cyls[best_cyl*7+2];
                float rad = rx*Dx + ry*Dy + rz*Dz;
                float n[3] = { rx - rad*Dx, ry - rad*Dy, rz - rad*Dz };
                float nn = rsqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] + 1e-30f);
                n[0]*=nn; n[1]*=nn; n[2]*=nn;
                reflectDir(dir, n);
                for (int i = 0; i < 3; ++i) { pr[i] += EPS*dir[i]; pv[i] += EPS*dir[i]; }
            } else if (best_axis >= 0) {                 // voxel face
                if (periodic) {
                    pv[best_axis] = (dir[best_axis] > 0.0f) ? vmin[best_axis] : vmax[best_axis];
                } else {
                    dir[best_axis] = -dir[best_axis];
                    pr[best_axis] += EPS*dir[best_axis]; pv[best_axis] += EPS*dir[best_axis];
                }
            } else {
                remaining = 0.0f;
            }
        }
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
    }
}

//! Moller-Trumbore ray-triangle intersection (port of Triangle::rayIntersects_MT).
//! Returns true and the signed distance t along (unit) dir to triangle (a,b,c).
__device__ inline bool mtHit(const float O[3], const float d[3],
                             const float a[3], const float b[3], const float c[3], float& t)
{
    const float DET_EPS = 1e-12f;
    const float TE = 1e-5f;            // barycentric tolerance (float32; generous to avoid edge leaks)
    float e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    float e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    float pv[3] = { d[1]*e2[2]-d[2]*e2[1], d[2]*e2[0]-d[0]*e2[2], d[0]*e2[1]-d[1]*e2[0] };
    float det = e1[0]*pv[0] + e1[1]*pv[1] + e1[2]*pv[2];
    if (det > -DET_EPS && det < DET_EPS) return false;
    float inv = 1.0f/det;
    float tv[3] = { O[0]-a[0], O[1]-a[1], O[2]-a[2] };
    float u = (tv[0]*pv[0] + tv[1]*pv[1] + tv[2]*pv[2]) * inv;
    if (u < -TE || u > 1.0f+TE) return false;
    float qv[3] = { tv[1]*e1[2]-tv[2]*e1[1], tv[2]*e1[0]-tv[0]*e1[2], tv[0]*e1[1]-tv[1]*e1[0] };
    float v = (d[0]*qv[0] + d[1]*qv[1] + d[2]*qv[2]) * inv;
    if (v < -TE || u+v > 1.0f+TE) return false;
    t = (e2[0]*qv[0] + e2[1]*qv[1] + e2[2]*qv[2]) * inv;
    return true;
}

//! Phase 2d/2e: one thread per walker, restricted diffusion inside ONE triangle mesh
//! within a voxel. Two positions (pv wrapped / pr unwrapped); each step takes the
//! NEAREST of {triangle crossing, voxel face}. Triangles reflect; voxel faces wrap pv
//! (periodic) or reflect (mirror). Walkers start uniformly in the voxel box (placement
//! intra/extra in a mesh would need a point-in-mesh test -- not done yet). No grid.
__global__ void meshKernel(float* traj, unsigned N, unsigned T, float step_length,
                           const float* verts, const int* tris, unsigned n_tris,
                           float vminx, float vminy, float vminz,
                           float vmaxx, float vmaxy, float vmaxz,
                           int periodic, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    float vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    float pv[3] = { vmin[0] + curand_uniform(&st)*(vmax[0]-vmin[0]),
                    vmin[1] + curand_uniform(&st)*(vmax[1]-vmin[1]),
                    vmin[2] + curand_uniform(&st)*(vmax[2]-vmin[2]) };
    float pr[3] = { pv[0], pv[1], pv[2] };

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];

    const float EPS = 1e-5f * step_length;   // re-hit guard / nudge (relative to step)
    for (unsigned t = 1; t <= T; ++t) {
        float gx = curand_normal(&st), gy = curand_normal(&st), gz = curand_normal(&st);
        float inv = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-30f);
        float dir[3] = { gx*inv, gy*inv, gz*inv };
        float remaining = step_length;

        for (int b = 0; b < 24 && remaining > 0.0f; ++b) {
            float O[3] = { pv[0], pv[1], pv[2] };
            float best_t = remaining; int best_tri = -1, best_axis = -1;
            for (unsigned i = 0; i < n_tris; ++i) {     // triangles (in pv)
                int i0 = tris[3*i+0], i1 = tris[3*i+1], i2 = tris[3*i+2];
                float a[3] = { verts[3*i0+0], verts[3*i0+1], verts[3*i0+2] };
                float bb[3]= { verts[3*i1+0], verts[3*i1+1], verts[3*i1+2] };
                float cc[3]= { verts[3*i2+0], verts[3*i2+1], verts[3*i2+2] };
                float tc;
                if (mtHit(O, dir, a, bb, cc, tc) && tc > EPS && tc < best_t) { best_t = tc; best_tri = int(i); best_axis = -1; }
            }
            for (int i = 0; i < 3; ++i) {               // voxel faces (in pv)
                float tw = -1.0f;
                if (dir[i] >  1e-20f) tw = (vmax[i]-pv[i])/dir[i];
                else if (dir[i] < -1e-20f) tw = (vmin[i]-pv[i])/dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; best_axis = i; best_tri = -1; }
            }
            for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]; pv[i] += best_t*dir[i]; }
            remaining -= best_t;
            if (best_tri >= 0) {                         // reflect off the triangle
                int i0 = tris[3*best_tri+0], i1 = tris[3*best_tri+1], i2 = tris[3*best_tri+2];
                float e1[3] = { verts[3*i1+0]-verts[3*i0+0], verts[3*i1+1]-verts[3*i0+1], verts[3*i1+2]-verts[3*i0+2] };
                float e2[3] = { verts[3*i2+0]-verts[3*i0+0], verts[3*i2+1]-verts[3*i0+1], verts[3*i2+2]-verts[3*i0+2] };
                float n[3] = { e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0] };
                float nn = rsqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] + 1e-30f);
                n[0]*=nn; n[1]*=nn; n[2]*=nn;
                reflectDir(dir, n);
                for (int i = 0; i < 3; ++i) { pr[i] += EPS*dir[i]; pv[i] += EPS*dir[i]; }
            } else if (best_axis >= 0) {                 // voxel face
                if (periodic) {
                    pv[best_axis] = (dir[best_axis] > 0.0f) ? vmin[best_axis] : vmax[best_axis];
                } else {
                    dir[best_axis] = -dir[best_axis];
                    pr[best_axis] += EPS*dir[best_axis]; pv[best_axis] += EPS*dir[best_axis];
                }
            } else {
                remaining = 0.0f;
            }
        }
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
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
    freeDiffusionKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length,
                                    p.vmin[0], p.vmin[1], p.vmin[2],
                                    p.vmax[0], p.vmax[1], p.vmax[2],
                                    p.periodic, p.seed);
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
    spheresKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, d_sph, n_sph,
                              p.vmin[0], p.vmin[1], p.vmin[2], p.vmax[0], p.vmax[1], p.vmax[2],
                              p.periodic, p.placement, p.seed);
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
    cylindersKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, d_cyl, n_cyl,
                                p.vmin[0], p.vmin[1], p.vmin[2], p.vmax[0], p.vmax[1], p.vmax[2],
                                p.periodic, p.placement, p.seed);
    CUDA_CHECK(cudaGetLastError());

    int rc = gpuSynthAndFree(p, d_traj, scheme, out_DWI);
    cudaFree(d_traj);
    cudaFree(d_cyl);
    return rc;
}

int runGpuMesh(const GpuSimParams& p, const float* verts, unsigned n_verts,
               const int* tris, unsigned n_tris,
               const float vmin[3], const float vmax[3],
               const float* scheme, double* out_DWI)
{
    float *d_traj = nullptr, *d_verts = nullptr;
    int   *d_tris = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj,  (size_t)p.N * (p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_verts, (size_t)n_verts * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_tris,  (size_t)n_tris  * 3 * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_verts, verts, (size_t)n_verts * 3 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tris,  tris,  (size_t)n_tris  * 3 * sizeof(int),   cudaMemcpyHostToDevice));

    const int bs = 256, gs = int((p.N + bs - 1) / bs);
    meshKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length, d_verts, d_tris, n_tris,
                           vmin[0], vmin[1], vmin[2], vmax[0], vmax[1], vmax[2], p.periodic, p.seed);
    CUDA_CHECK(cudaGetLastError());

    int rc = gpuSynthAndFree(p, d_traj, scheme, out_DWI);
    cudaFree(d_traj);
    cudaFree(d_verts);
    cudaFree(d_tris);
    return rc;
}
