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
#include <cstdlib>
#include <chrono>
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
                                    float sminx, float sminy, float sminz,
                                    float smaxx, float smaxy, float smaxz,
                                    int periodic, float dt, float t2_extra,
                                    int* out_inobs, float* t2sum, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);
    out_inobs[w] = -1;                              // free diffusion: every walker is extra

    float vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    float smin[3] = { sminx, sminy, sminz }, smax[3] = { smaxx, smaxy, smaxz };
    // start uniformly in the sampling area; pr and pv coincide at t=0
    float pv[3] = { smin[0] + curand_uniform(&st)*(smax[0]-smin[0]),
                    smin[1] + curand_uniform(&st)*(smax[1]-smin[1]),
                    smin[2] + curand_uniform(&st)*(smax[2]-smin[2]) };
    float pr[3] = { pv[0], pv[1], pv[2] };
    // Per-axis mirror sign: pos_v moves by dir, pos_r by dir*bmir. Mirror face flips both ->
    // pos_r continues straight (free, correct MSD); pos_v reflects in. Periodic keeps bmir=+1.
    float bmir[3] = { 1.0f, 1.0f, 1.0f };
    float voxext  = fmaxf(fmaxf(vmax[0]-vmin[0], vmax[1]-vmin[1]), vmax[2]-vmin[2]);

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];
    size_t tb = (size_t)w * (T + 1); t2sum[tb] = 0.0f; float t2acc = 0.0f;  // Phase 6b (free: all extra)

    const float NUDGE = fmaxf(1e-5f * step_length, 1e-5f * voxext);   // float32-safe wall nudge
    const float EPS   = 0.25f * NUDGE;
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
                for (int i = 0; i < 3; ++i) { pr[i] += remaining*dir[i]*bmir[i]; pv[i] += remaining*dir[i]; }
                remaining = 0.0f;
            } else {
                for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]*bmir[i]; pv[i] += best_t*dir[i]; }
                remaining -= best_t;
                if (periodic) {
                    pv[axis] = (dir[axis] > 0.0f) ? vmin[axis] : vmax[axis];   // wrap pv only; pr continuous
                } else {
                    // MIRROR: flip dir AND bmir -> pv reflects in, pr keeps its real direction (free).
                    dir[axis]  = -dir[axis];
                    bmir[axis] = -bmir[axis];
                    pv[axis] += NUDGE*dir[axis];
                    pr[axis] += NUDGE*dir[axis]*bmir[axis];
                }
            }
        }
        t2acc += dt / t2_extra; t2sum[tb + t] = t2acc;     // Phase 6b: free walker is always extra
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
// Index of the sphere containing (x,y,z), or -1 if outside all of them (= extra).
__device__ inline int whichSphere(float x, float y, float z, const float* spheres, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        float dx = x - spheres[i*4+0], dy = y - spheres[i*4+1], dz = z - spheres[i*4+2];
        float r = spheres[i*4+3];
        if (dx*dx + dy*dy + dz*dz < r*r) return int(i);
    }
    return -1;
}
__device__ inline bool insideAnySphere(float x, float y, float z, const float* spheres, unsigned n)
{
    return whichSphere(x, y, z, spheres, n) >= 0;
}

//! Phase 2a/2c/2e: restricted diffusion among impermeable spheres inside a voxel.
//! Two positions: pv (wrapped, for collision) and pr (unwrapped, for the signal).
//! Each step takes the nearest of {sphere surfaces, voxel faces}; spheres reflect,
//! voxel faces wrap pv (periodic) or reflect (mirror, folded). placement: 0=uniform
//! in voxel, 1=intra (random sphere), 2=extra (voxel, outside all spheres).
__global__ void spheresKernel(float* traj, unsigned N, unsigned T, float step_intra, float step_extra,
                              const float* spheres, unsigned n_sph,
                              const float* prob_ie, const float* prob_ei,
                              const int* cell_start, const int* cell_sph,   // CSR spatial grid (per-cell sphere lists)
                              double gminx, double gminy, double gminz,     // grid origin
                              double icsx, double icsy, double icsz,        // inverse cell size (cells per unit)
                              int gdimx, int gdimy, int gdimz,              // grid dims
                              float vminx, float vminy, float vminz,
                              float vmaxx, float vmaxy, float vmaxz,
                              float sminx, float sminy, float sminz,
                              float smaxx, float smaxy, float smaxz,
                              int periodic, int placement, float dt, float t2_extra,
                              const float* obs_T2, int* out_inobs, float* t2sum,
                              unsigned long long seed, unsigned w_offset)
{
    const double gmin[3] = { gminx, gminy, gminz }, ics[3] = { icsx, icsy, icsz };
    const int    gdim[3] = { gdimx, gdimy, gdimz };
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;       // N is the CHUNK size; traj/t2sum/out_inobs are chunk-local

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w_offset + w, 0, &st);   // global walker id -> chunk-size-independent

    // Geometry in DOUBLE: float32 cannot place a walker exactly on a curved membrane, so a
    // float reflect lands a hair off-surface and the next step tunnels through -> a per-bounce
    // leak that accumulates over many bounces (tight packings / long times). Double matches the
    // CPU exactly, so the reflect is clean and NO nudge / containment hack is needed.
    double vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    double smin[3] = { sminx, sminy, sminz }, smax[3] = { smaxx, smaxy, smaxz };
    double pv[3]; int in_obs;
    if (placement == 1) {                          // intra: uniform in a random sphere
        unsigned si = curand(&st) % n_sph;
        double cx = spheres[si*4+0], cy = spheres[si*4+1], cz = spheres[si*4+2], r0 = spheres[si*4+3];
        double ox, oy, oz;
        do { ox = (2.0*curand_uniform(&st)-1.0)*r0; oy = (2.0*curand_uniform(&st)-1.0)*r0; oz = (2.0*curand_uniform(&st)-1.0)*r0; }
        while (ox*ox + oy*oy + oz*oz > r0*r0);
        pv[0] = cx + ox; pv[1] = cy + oy; pv[2] = cz + oz; in_obs = int(si);
    } else {                                       // uniform (0) or extra (2): in the sampling area
        for (int tries = 0; tries < 10000; ++tries) {
            pv[0] = smin[0] + (double)curand_uniform(&st)*(smax[0]-smin[0]);
            pv[1] = smin[1] + (double)curand_uniform(&st)*(smax[1]-smin[1]);
            pv[2] = smin[2] + (double)curand_uniform(&st)*(smax[2]-smin[2]);
            if (placement != 2 || !insideAnySphere(pv[0], pv[1], pv[2], spheres, n_sph)) break;
        }
        in_obs = whichSphere(pv[0], pv[1], pv[2], spheres, n_sph);   // compartment at placement
    }
    double pr[3] = { pv[0], pv[1], pv[2] };
    // Per-axis mirror sign (CPU boundary_mirror): pos_v moves by dir, pos_r by dir*bmir.
    // Mirror face flips both -> pos_r continuous (free MSD); pos_v reflects in. Periodic: bmir=+1.
    double bmir[3] = { 1.0, 1.0, 1.0 };
    double voxext  = fmax(fmax(vmax[0]-vmin[0], vmax[1]-vmin[1]), vmax[2]-vmin[2]);

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];
    size_t tb = (size_t)w * (T + 1); t2sum[tb] = 0.0f; float t2acc = 0.0f;  // Phase 6b

    for (unsigned t = 1; t <= T; ++t) {
        // per-timestep step length from the CURRENT compartment (permeable crossings may
        // have changed in_obs). EPS just skips the ~0 root of the surface we rest on.
        double step_length = (in_obs >= 0) ? step_intra : step_extra;
        double EPS = 1e-9 * voxext;                      // << any real gap; double-clean

        double dx = curand_normal(&st), dy = curand_normal(&st), dz = curand_normal(&st);
        double inv = 1.0 / sqrt(dx*dx + dy*dy + dz*dz + 1e-300);
        double dir[3] = { dx*inv, dy*inv, dz*inv };
        double remaining = step_length;
        int rest_sph = -1;                               // sphere we are currently resting on

        for (int b = 0; b < 24 && remaining > EPS; ++b) {
            double best_t = remaining; int best_sph = -1, best_axis = -1;
            // Only test spheres in the grid cells the step segment [pv, pv+remaining*dir] overlaps
            // (clamped). A real hit's contact point lies in BOTH the segment box and the sphere's
            // AABB, so it shares a queried cell -> no missed collision (same guarantee as the mesh
            // grid). O(nearby) instead of O(n_sph); for tiny lists the grid is ~1 cell = brute.
            double segE[3] = { pv[0]+remaining*dir[0], pv[1]+remaining*dir[1], pv[2]+remaining*dir[2] };
            int cmn[3], cmx[3];
            for (int a = 0; a < 3; ++a) {
                double lo3 = fmin(pv[a], segE[a]) - gmin[a], hi3 = fmax(pv[a], segE[a]) - gmin[a];
                int clo = int(floor(lo3*ics[a])), chi = int(floor(hi3*ics[a]));
                cmn[a] = clo < 0 ? 0 : (clo > gdim[a]-1 ? gdim[a]-1 : clo);
                cmx[a] = chi < 0 ? 0 : (chi > gdim[a]-1 ? gdim[a]-1 : chi);
            }
            for (int cz = cmn[2]; cz <= cmx[2]; ++cz)
            for (int cy = cmn[1]; cy <= cmx[1]; ++cy)
            for (int cx = cmn[0]; cx <= cmx[0]; ++cx) {
                int cell = cx + gdim[0]*(cy + gdim[1]*cz);
                for (int kk = cell_start[cell]; kk < cell_start[cell+1]; ++kk) {
                    unsigned i = (unsigned)cell_sph[kk];     // sphere surfaces (collide in pv)
                    double ux = pv[0]-(double)spheres[i*4+0], uy = pv[1]-(double)spheres[i*4+1], uz = pv[2]-(double)spheres[i*4+2];
                    double ri = spheres[i*4+3];
                    double bq = ux*dir[0] + uy*dir[1] + uz*dir[2];
                    double cq = ux*ux + uy*uy + uz*uz - ri*ri;
                    double discr = bq*bq - cq;
                    if (discr > 0.0) {
                        double sq = sqrt(discr);
                        double t1 = -bq - sq, t2 = -bq + sq;
                        // Honor the nearest FORWARD root for every sphere we are NOT resting on (never
                        // skip a real approaching entry -- a blanket tc>EPS skip tunnels to the far root
                        // and leaks). For the surface we rest on, skip its ~0 root so we don't re-hit it.
                        double lo = (int(i) == rest_sph) ? EPS : 0.0;
                        double tc = (t1 > lo) ? t1 : ((t2 > lo) ? t2 : -1.0);
                        if (tc > lo && tc < best_t) { best_t = tc; best_sph = int(i); best_axis = -1; }
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {                // voxel faces (in pv)
                double tw = -1.0;
                if (dir[i] >  1e-30) tw = (vmax[i]-pv[i])/dir[i];
                else if (dir[i] < -1e-30) tw = (vmin[i]-pv[i])/dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; best_axis = i; best_sph = -1; }
            }
            for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]*bmir[i]; pv[i] += best_t*dir[i]; }
            remaining -= best_t;
            if (best_sph >= 0) {                          // sphere membrane: cross or reflect
                bool from_intra = (in_obs == best_sph);
                float prob = from_intra ? prob_ie[best_sph] : prob_ei[best_sph];
                if (prob > 0.0f && prob > curand_uniform(&st)) {
                    // PERMEATE: transparent membrane -- finish straight through, re-classify.
                    for (int i = 0; i < 3; ++i) { pr[i] += remaining*dir[i]*bmir[i]; pv[i] += remaining*dir[i]; }
                    remaining = 0.0;
                    in_obs = whichSphere(pv[0], pv[1], pv[2], spheres, n_sph);
                } else {                                  // reflect off the sphere (double => lands exactly on the surface; no nudge)
                    double nx = pv[0]-(double)spheres[best_sph*4+0], ny = pv[1]-(double)spheres[best_sph*4+1], nz = pv[2]-(double)spheres[best_sph*4+2];
                    double nn = 1.0 / sqrt(nx*nx + ny*ny + nz*nz + 1e-300);
                    double n[3] = { nx*nn, ny*nn, nz*nn };
                    double dn = dir[0]*n[0] + dir[1]*n[1] + dir[2]*n[2];
                    dir[0] -= 2.0*dn*n[0]; dir[1] -= 2.0*dn*n[1]; dir[2] -= 2.0*dn*n[2];
                    rest_sph = best_sph;                  // now resting on this surface (its ~0 root is skipped)
                }
            } else if (best_axis >= 0) {                 // voxel face
                if (periodic) {
                    pv[best_axis] = (dir[best_axis] > 0.0) ? vmin[best_axis] : vmax[best_axis];
                } else {
                    dir[best_axis]  = -dir[best_axis];
                    bmir[best_axis] = -bmir[best_axis];
                }
                rest_sph = -1;                            // crossed a wall: no resting sphere
            } else {
                remaining = 0.0;                          // no event: full step taken
            }
        }
        in_obs = whichSphere(pv[0], pv[1], pv[2], spheres, n_sph);  // compartment for next step
        float T2c = (in_obs >= 0) ? obs_T2[in_obs] : t2_extra;
        t2acc += dt / T2c; t2sum[tb + t] = t2acc;
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
    }
    // Re-classify the landing position (matches CPU finalPositionCheck, which relabels by
    // walker.location = isInIntra(pos_v)). Analytic obstacles are NOT deported on the CPU (its
    // deport guard requires plyObstacles_list), so we relabel but do NOT discard (-2), keeping
    // GPU and CPU consistent. Relabel leaves the total signal unchanged, only the intra/extra split.
    out_inobs[w] = whichSphere(pv[0], pv[1], pv[2], spheres, n_sph);
}

//! Phase 2b/2c: one thread per walker, restricted diffusion among impermeable
//! infinite cylinders (cyls[i*7+0..6] = px,py,pz, dx,dy,dz(unit), r). Walker starts
//! in a random cylinder's disk; each step takes the NEAREST radial crossing over
//! all cylinders and reflects off that radial normal (axial motion is free).
// Index of the cylinder containing (x,y,z) (radially), or -1 if outside all (= extra).
__device__ inline int whichCylinder(float x, float y, float z, const float* cyls, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        float ex = x-cyls[i*7+0], ey = y-cyls[i*7+1], ez = z-cyls[i*7+2];
        float Dx = cyls[i*7+3], Dy = cyls[i*7+4], Dz = cyls[i*7+5], r = cyls[i*7+6];
        float md = ex*Dx + ey*Dy + ez*Dz;
        float rad2 = ex*ex + ey*ey + ez*ez - md*md;
        if (rad2 < r*r) return int(i);
    }
    return -1;
}
__device__ inline bool insideAnyCylinder(float x, float y, float z, const float* cyls, unsigned n)
{
    return whichCylinder(x, y, z, cyls, n) >= 0;
}

//! Phase 2b/2c/2e: restricted diffusion among impermeable infinite cylinders inside
//! a voxel. Two positions (pv wrapped / pr unwrapped); each step takes the nearest of
//! {cylinder surfaces, voxel faces}. placement: 0=uniform, 1=intra (random cylinder
//! disk), 2=extra (in the voxel, outside all cylinders).
__global__ void cylindersKernel(float* traj, unsigned N, unsigned T, float step_intra, float step_extra,
                                const float* cyls, unsigned n_cyl,
                                const float* prob_ie, const float* prob_ei,
                                const int* cell_start, const int* cell_cyl,   // CSR spatial grid
                                double gminx, double gminy, double gminz,
                                double icsx, double icsy, double icsz,
                                int gdimx, int gdimy, int gdimz,
                                float vminx, float vminy, float vminz,
                                float vmaxx, float vmaxy, float vmaxz,
                                float sminx, float sminy, float sminz,
                                float smaxx, float smaxy, float smaxz,
                                int periodic, int placement, float dt, float t2_extra,
                                const float* obs_T2, int* out_inobs, float* t2sum, unsigned long long seed)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;
    const double gmin[3] = { gminx, gminy, gminz }, ics[3] = { icsx, icsy, icsz };
    const int    gdim[3] = { gdimx, gdimy, gdimz };

    curandStatePhilox4_32_10_t st;
    curand_init(seed, w, 0, &st);

    // Geometry in DOUBLE (see spheresKernel): float32 cannot place a walker exactly on a curved
    // membrane -> per-bounce micro-leak that accumulates over many bounces. Double matches the CPU.
    double vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    double smin[3] = { sminx, sminy, sminz }, smax[3] = { smaxx, smaxy, smaxz };
    double pv[3]; int in_obs;
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
        in_obs = int(ci);
    } else {                                       // uniform (0) or extra (2): in the sampling area
        for (int tries = 0; tries < 10000; ++tries) {
            pv[0] = smin[0] + curand_uniform(&st)*(smax[0]-smin[0]);
            pv[1] = smin[1] + curand_uniform(&st)*(smax[1]-smin[1]);
            pv[2] = smin[2] + curand_uniform(&st)*(smax[2]-smin[2]);
            if (placement != 2 || !insideAnyCylinder(pv[0], pv[1], pv[2], cyls, n_cyl)) break;
        }
        in_obs = whichCylinder(pv[0], pv[1], pv[2], cyls, n_cyl);
    }
    double pr[3] = { pv[0], pv[1], pv[2] };
    // Per-axis mirror sign (CPU boundary_mirror): pos_v by dir, pos_r by dir*bmir. Mirror face
    // flips both -> pos_r continuous (free MSD); pos_v reflects in. Periodic keeps bmir=+1.
    double bmir[3] = { 1.0, 1.0, 1.0 };
    double voxext  = fmax(fmax(vmax[0]-vmin[0], vmax[1]-vmin[1]), vmax[2]-vmin[2]);

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];
    size_t tb = (size_t)w * (T + 1); t2sum[tb] = 0.0f; float t2acc = 0.0f;  // Phase 6b

    for (unsigned t = 1; t <= T; ++t) {
        double step_length = (in_obs >= 0) ? step_intra : step_extra;  // current compartment (permeability)
        double EPS = 1e-9 * voxext;                       // skips the ~0 root of the surface we rest on
        double gx = curand_normal(&st), gy = curand_normal(&st), gz = curand_normal(&st);
        double inv = 1.0 / sqrt(gx*gx + gy*gy + gz*gz + 1e-300);
        double dir[3] = { gx*inv, gy*inv, gz*inv };
        double remaining = step_length;
        int rest_cyl = -1;                                // cylinder we are currently resting on

        for (int b = 0; b < 24 && remaining > EPS; ++b) {
            double best_t = remaining; int best_cyl = -1, best_axis = -1;
            // Only test cylinders in the grid cells the step segment overlaps (clamped) -- same
            // O(nearby) query + no-missed-collision guarantee as the sphere grid. Axis-aligned
            // (z-parallel) packings prune well; tilted cylinders span the grid (=> ~brute, correct).
            double segE[3] = { pv[0]+remaining*dir[0], pv[1]+remaining*dir[1], pv[2]+remaining*dir[2] };
            int cmn[3], cmx[3];
            for (int aa = 0; aa < 3; ++aa) {
                double lo3 = fmin(pv[aa], segE[aa]) - gmin[aa], hi3 = fmax(pv[aa], segE[aa]) - gmin[aa];
                int clo = int(floor(lo3*ics[aa])), chi = int(floor(hi3*ics[aa]));
                cmn[aa] = clo < 0 ? 0 : (clo > gdim[aa]-1 ? gdim[aa]-1 : clo);
                cmx[aa] = chi < 0 ? 0 : (chi > gdim[aa]-1 ? gdim[aa]-1 : chi);
            }
            for (int cz = cmn[2]; cz <= cmx[2]; ++cz)
            for (int cy = cmn[1]; cy <= cmx[1]; ++cy)
            for (int cx = cmn[0]; cx <= cmx[0]; ++cx) {
                int cell = cx + gdim[0]*(cy + gdim[1]*cz);
                for (int kk = cell_start[cell]; kk < cell_start[cell+1]; ++kk) {
                    unsigned i = (unsigned)cell_cyl[kk];   // cylinder surfaces (in pv)
                    double ex = pv[0]-(double)cyls[i*7+0], ey = pv[1]-(double)cyls[i*7+1], ez = pv[2]-(double)cyls[i*7+2];
                    double Dx = cyls[i*7+3], Dy = cyls[i*7+4], Dz = cyls[i*7+5], ri = cyls[i*7+6];
                    double md = ex*Dx + ey*Dy + ez*Dz;
                    double nd = dir[0]*Dx + dir[1]*Dy + dir[2]*Dz;
                    double a  = 1.0 - nd*nd;
                    if (a <= 1e-12) continue;
                    double cq = (ex*ex + ey*ey + ez*ez - ri*ri) - md*md;
                    double bq = (ex*dir[0] + ey*dir[1] + ez*dir[2]) - nd*md;
                    double discr = bq*bq - a*cq;
                    if (discr > 0.0) {
                        double sq = sqrt(discr);
                        double t1 = (-bq - sq) / a, t2 = (-bq + sq) / a;
                        // Honor the nearest forward root for cylinders we are NOT resting on; skip the
                        // ~0 root only for the one we rest on (no blanket EPS skip -> no tunnel).
                        double lo = (int(i) == rest_cyl) ? EPS : 0.0;
                        double tc = (t1 > lo) ? t1 : ((t2 > lo) ? t2 : -1.0);
                        if (tc > lo && tc < best_t) { best_t = tc; best_cyl = int(i); best_axis = -1; }
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {               // voxel faces (in pv)
                double tw = -1.0;
                if (dir[i] >  1e-30) tw = (vmax[i]-pv[i])/dir[i];
                else if (dir[i] < -1e-30) tw = (vmin[i]-pv[i])/dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; best_axis = i; best_cyl = -1; }
            }
            for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]*bmir[i]; pv[i] += best_t*dir[i]; }
            remaining -= best_t;
            if (best_cyl >= 0) {                         // cylinder membrane: cross or reflect
                bool from_intra = (in_obs == best_cyl);
                float prob = from_intra ? prob_ie[best_cyl] : prob_ei[best_cyl];
                if (prob > 0.0f && prob > curand_uniform(&st)) {
                    // PERMEATE (CPU model): finish the step straight through, re-classify at landing.
                    for (int i = 0; i < 3; ++i) { pr[i] += remaining*dir[i]*bmir[i]; pv[i] += remaining*dir[i]; }
                    remaining = 0.0;
                    in_obs = whichCylinder(pv[0], pv[1], pv[2], cyls, n_cyl);
                } else {                                  // reflect off the radial normal (double => exact; no nudge)
                    double Dx = cyls[best_cyl*7+3], Dy = cyls[best_cyl*7+4], Dz = cyls[best_cyl*7+5];
                    double rx = pv[0]-(double)cyls[best_cyl*7+0], ry = pv[1]-(double)cyls[best_cyl*7+1], rz = pv[2]-(double)cyls[best_cyl*7+2];
                    double rad = rx*Dx + ry*Dy + rz*Dz;
                    double n[3] = { rx - rad*Dx, ry - rad*Dy, rz - rad*Dz };
                    double nn = 1.0 / sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] + 1e-300);
                    n[0]*=nn; n[1]*=nn; n[2]*=nn;
                    double dn = dir[0]*n[0] + dir[1]*n[1] + dir[2]*n[2];
                    dir[0] -= 2.0*dn*n[0]; dir[1] -= 2.0*dn*n[1]; dir[2] -= 2.0*dn*n[2];
                    rest_cyl = best_cyl;
                }
            } else if (best_axis >= 0) {                 // voxel face
                if (periodic) {
                    pv[best_axis] = (dir[best_axis] > 0.0) ? vmin[best_axis] : vmax[best_axis];   // pr continuous
                } else {
                    // MIRROR: flip dir AND bmir -> pv reflects in, pr keeps real direction (free).
                    dir[best_axis]  = -dir[best_axis];
                    bmir[best_axis] = -bmir[best_axis];
                }
                rest_cyl = -1;                            // crossed a wall: no resting cylinder
            } else {
                remaining = 0.0;
            }
        }
        float T2c = (in_obs >= 0) ? obs_T2[in_obs] : t2_extra;   // T2 of the END-of-step compartment
        t2acc += dt / T2c; t2sum[tb + t] = t2acc;                // Phase 6b: per-step T2 accumulation
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
    }
    // Re-classify the landing position (matches CPU finalPositionCheck). Analytic obstacles are
    // NOT deported on the CPU, so relabel but do NOT discard -- keeps GPU/CPU consistent.
    out_inobs[w] = whichCylinder(pv[0], pv[1], pv[2], cyls, n_cyl);
}

//! Moller-Trumbore ray-triangle intersection (port of Triangle::rayIntersects_MT).
//! Returns true and the signed distance t along (unit) dir to triangle (a,b,c).
// Möller–Trumbore in DOUBLE: float32 cannot place a walker exactly on a triangle membrane,
// so a float reflect lands a hair off-surface and the next step tunnels through -> a per-bounce
// leak that accumulates on thin / small-triangle meshes (e.g. a 1um-radius tube). Double matches
// the CPU exactly, so the reflect is clean and NO nudge / containment hack is needed.
__device__ inline bool mtHit(const double O[3], const double d[3],
                             const double a[3], const double b[3], const double c[3], double& t)
{
    const double DET_EPS = 1e-24;
    const double TE = 1e-7;            // barycentric tolerance (double; small sliver, still edge-leak-safe)
    double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    double e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    double pv[3] = { d[1]*e2[2]-d[2]*e2[1], d[2]*e2[0]-d[0]*e2[2], d[0]*e2[1]-d[1]*e2[0] };
    double det = e1[0]*pv[0] + e1[1]*pv[1] + e1[2]*pv[2];
    if (det > -DET_EPS && det < DET_EPS) return false;
    double inv = 1.0/det;
    double tv[3] = { O[0]-a[0], O[1]-a[1], O[2]-a[2] };
    double u = (tv[0]*pv[0] + tv[1]*pv[1] + tv[2]*pv[2]) * inv;
    if (u < -TE || u > 1.0+TE) return false;
    double qv[3] = { tv[1]*e1[2]-tv[2]*e1[1], tv[2]*e1[0]-tv[0]*e1[2], tv[0]*e1[1]-tv[1]*e1[0] };
    double v = (d[0]*qv[0] + d[1]*qv[1] + d[2]*qv[2]) * inv;
    if (v < -TE || u+v > 1.0+TE) return false;
    t = (e2[0]*qv[0] + e2[1]*qv[1] + e2[2]*qv[2]) * inv;
    return true;
}

//! Phase 3b: classify a point against the merged mesh -> the owning obstacle id, or
//! -1 (extra), by RAY CASTING (parity), the robust point-in-mesh test. A nearest-triangle
//! test is WRONG: a large triangle whose centroid is far can still be the containing
//! surface while a small nearby triangle's normal is unrepresentative.
//!
//! We shoot a +x ray from P. Because it keeps y,z fixed, every crossing lies in the same
//! grid cell-row (cy,cz)=P's cell, so we only march that row (cheap and exact). Each
//! triangle is registered in every cell its AABB overlaps; to count a crossing once we
//! accept it only in the cell that contains the crossing's x (standard grid dedup).
//!  - No forward crossing  -> P is outside everything -> extra (correct for sparse meshes
//!    / points beyond the mesh; degrades gracefully instead of guessing a far triangle).
//!  - Otherwise P is inside obstacle m == owner of the NEAREST crossing IFF the ray crosses
//!    m's surface an ODD number of times (parity handles concave shells). Obstacles do not
//!    overlap, so at most one has odd parity, and it is the nearest crossing's owner.
//! Two short marches of the row (find nearest; count that owner's parity) keep it
//! register-only (no per-obstacle scratch), so it scales to many obstacles.
__device__ inline int classifyMesh(const double P[3], const float* verts, const int* tris,
                                   const int* tri_obs,
                                   const double gmin[3], const double ics[3], const int dim[3],
                                   const int* cell_start, const int* cell_tri)
{
    int pcx = int(floor((P[0]-gmin[0])*ics[0]));    // start cell x (clamped); ray goes +x
    if (pcx < 0) pcx = 0; else if (pcx > dim[0]-1) pcx = dim[0]-1;
    int cy = int(floor((P[1]-gmin[1])*ics[1]));
    int cz = int(floor((P[2]-gmin[2])*ics[2]));
    if (cy < 0 || cy >= dim[1] || cz < 0 || cz >= dim[2]) return -1;   // ray row outside grid
    // Sub-nm forward-crossing cutoff. The old float32 EPS=1e-7 (=0.1um) skipped crossings
    // within 0.1um ahead, so a walker near the inner wall of a ~1um object had its wall
    // crossing dropped and was misclassified as extra (-> b0 loss + false leak flags). In
    // double, EPS just excludes the crossing at the point itself.
    const double EPS = 1e-12;
    const int   rowoff = dim[0]*(cy + dim[1]*cz);

    // March 1: nearest forward +x crossing (dedup: count only in the cell holding its x).
    double best_t = 1e30; int best_obs = -1;
    for (int cx = pcx; cx < dim[0]; ++cx) {
        int cell = cx + rowoff;
        for (int k = cell_start[cell]; k < cell_start[cell+1]; ++k) {
            int i = cell_tri[k];
            int i0 = tris[3*i+0], i1 = tris[3*i+1], i2 = tris[3*i+2];
            double a[3]  = { verts[3*i0+0], verts[3*i0+1], verts[3*i0+2] };
            double e1[3] = { verts[3*i1+0]-a[0], verts[3*i1+1]-a[1], verts[3*i1+2]-a[2] };
            double e2[3] = { verts[3*i2+0]-a[0], verts[3*i2+1]-a[1], verts[3*i2+2]-a[2] };
            // ray-triangle for d=(1,0,0): t = signed x-distance to the triangle plane/bary.
            double pvx = -e2[2], pvy = e2[1];                // d x e2 = (0, -e2z, e2y)
            double det = e1[1]*pvx + e1[2]*pvy;              // e1 . (d x e2)
            if (det > -1e-40 && det < 1e-40) continue;       // ray parallel to triangle
            double inv = 1.0/det;
            double tv[3] = { P[0]-a[0], P[1]-a[1], P[2]-a[2] };
            double u = (tv[1]*pvx + tv[2]*pvy) * inv;
            if (u < -1e-7 || u > 1.0+1e-7) continue;
            double qx = tv[1]*e1[2]-tv[2]*e1[1], qy = tv[2]*e1[0]-tv[0]*e1[2], qz = tv[0]*e1[1]-tv[1]*e1[0];
            double v = qx * inv;                             // d . q = qx
            if (v < -1e-7 || u+v > 1.0+1e-7) continue;
            double t = (e2[0]*qx + e2[1]*qy + e2[2]*qz) * inv;
            if (t <= EPS) continue;                          // backward / on the point
            int hitcx = int(floor((P[0]+t-gmin[0])*ics[0]));
            if (hitcx != cx) continue;                       // dedup: count in its own cell only
            if (t < best_t) { best_t = t; best_obs = tri_obs[i]; }
        }
        if (best_obs >= 0) break;   // nearest is in the lowest cx that has any crossing
    }
    if (best_obs < 0) return -1;    // no surface ahead -> extra

    // March 2: parity of best_obs's crossings along the whole +x row.
    int crossings = 0;
    for (int cx = pcx; cx < dim[0]; ++cx) {
        int cell = cx + rowoff;
        for (int k = cell_start[cell]; k < cell_start[cell+1]; ++k) {
            int i = cell_tri[k];
            if (tri_obs[i] != best_obs) continue;
            int i0 = tris[3*i+0], i1 = tris[3*i+1], i2 = tris[3*i+2];
            double a[3]  = { verts[3*i0+0], verts[3*i0+1], verts[3*i0+2] };
            double e1[3] = { verts[3*i1+0]-a[0], verts[3*i1+1]-a[1], verts[3*i1+2]-a[2] };
            double e2[3] = { verts[3*i2+0]-a[0], verts[3*i2+1]-a[1], verts[3*i2+2]-a[2] };
            double pvx = -e2[2], pvy = e2[1];
            double det = e1[1]*pvx + e1[2]*pvy;
            if (det > -1e-40 && det < 1e-40) continue;
            double inv = 1.0/det;
            double tv[3] = { P[0]-a[0], P[1]-a[1], P[2]-a[2] };
            double u = (tv[1]*pvx + tv[2]*pvy) * inv;
            if (u < -1e-7 || u > 1.0+1e-7) continue;
            double qx = tv[1]*e1[2]-tv[2]*e1[1], qy = tv[2]*e1[0]-tv[0]*e1[2], qz = tv[0]*e1[1]-tv[1]*e1[0];
            double v = qx * inv;
            if (v < -1e-7 || u+v > 1.0+1e-7) continue;
            double t = (e2[0]*qx + e2[1]*qy + e2[2]*qz) * inv;
            if (t <= EPS) continue;
            int hitcx = int(floor((P[0]+t-gmin[0])*ics[0]));
            if (hitcx != cx) continue;
            crossings++;
        }
    }
    return (crossings & 1) ? best_obs : -1;
}

//! Phase 2d/2e/2f: one thread per walker, restricted diffusion inside ONE triangle
//! mesh within a voxel, accelerated by a uniform grid (CSR cell->triangle lists).
//! Two positions (pv wrapped / pr unwrapped); each step takes the NEAREST of {triangle
//! crossing, voxel face}. Triangles reflect; voxel faces wrap pv (periodic) or reflect
//! (mirror). Instead of all triangles, only those in the grid cells the step segment
//! overlaps are tested (the hit cell is in the segment AABB, so this is exact).
//! Grid: origin gmin, inverse cell size ics, dims dim; cell c=(cx,cy,cz) holds the
//! triangle indices cell_tri[cell_start[c] .. cell_start[c+1]).
__global__ void meshKernel(float* traj, unsigned N, unsigned T, float step_extra,
                           const float* verts, const int* tris, unsigned n_tris,
                           const int* tri_obs, const float* obs_T2, const float* obs_step,
                           const float* prob_ie, const float* prob_ei,
                           float vminx, float vminy, float vminz,
                           float vmaxx, float vmaxy, float vmaxz, int periodic, int placement,
                           float sminx, float sminy, float sminz,
                           float smaxx, float smaxy, float smaxz,
                           double gminx, double gminy, double gminz,
                           double icsx, double icsy, double icsz,
                           int dimx, int dimy, int dimz,
                           const int* cell_start, const int* cell_tri,
                           float dt, float t2_extra, int* out_inobs, float* t2sum,
                           unsigned long long seed, unsigned w_offset)
{
    unsigned w = blockIdx.x * blockDim.x + threadIdx.x;
    if (w >= N) return;       // N is the CHUNK size; traj/t2sum/out_inobs are chunk-local (index w)

    curandStatePhilox4_32_10_t st;
    // Seed by the GLOBAL walker id (w_offset + w) so results are identical regardless of
    // how the N walkers are split into chunks (reproducible, chunk-size-independent).
    curand_init(seed, w_offset + w, 0, &st);

    double vmin[3] = { vminx, vminy, vminz }, vmax[3] = { vmaxx, vmaxy, vmaxz };
    double smin[3] = { sminx, sminy, sminz }, smax[3] = { smaxx, smaxy, smaxz };
    double gmin[3] = { gminx, gminy, gminz }, ics[3] = { icsx, icsy, icsz };
    int    dim[3]  = { dimx, dimy, dimz };
    // Geometry in DOUBLE (matches the CPU and the sphere/cylinder kernels): the reflect lands
    // exactly on the triangle, so NO off-surface nudge is needed. voxext scales the tiny EPS
    // that skips the ~0 root of the face the walker rests on.
    double voxext = fmax(fmax(vmax[0]-vmin[0], vmax[1]-vmin[1]), vmax[2]-vmin[2]);
    // Placement honoring ini_walkers_pos (CPU parity): rejection-sample in the sampling area
    // until the requested compartment is hit -- intra (1) = inside any mesh obstacle, extra
    // (2) = outside all; uniform (0) accepts the first sample. Without this the mesh kernel
    // seeded UNIFORMLY in the voxel regardless of the flag, so an "intra" run started only the
    // area-fraction of walkers inside the mesh (the rest free-diffused) -> the reported
    // analytic-vs-PLY mismatch. classifyMesh gives the membership (-1 = extra). in_obs then
    // selects the per-compartment step length / T2 and follows permeable crossings mid-walk.
    double pv[3]; int in_obs = -1;
    for (int tries = 0; tries < 10000; ++tries) {
        pv[0] = smin[0] + (double)curand_uniform(&st)*(smax[0]-smin[0]);
        pv[1] = smin[1] + (double)curand_uniform(&st)*(smax[1]-smin[1]);
        pv[2] = smin[2] + (double)curand_uniform(&st)*(smax[2]-smin[2]);
        in_obs = classifyMesh(pv, verts, tris, tri_obs, gmin, ics, dim, cell_start, cell_tri);
        if (placement == 0) break;                    // uniform: accept any position
        if (placement == 1 && in_obs >= 0) break;     // intra: inside a mesh obstacle
        if (placement == 2 && in_obs <  0) break;     // extra: outside all obstacles
    }
    double pr[3] = { pv[0], pv[1], pv[2] };
    // Per-axis mirror sign (CPU walker.boundary_mirror). pos_v moves by dir (voxel frame,
    // reflects/wraps at faces); pos_r moves by dir*bmir (real frame). A mirror face flips
    // BOTH dir[a] and bmir[a] -> their product is unchanged, so pos_r continues straight
    // (free, unbounded -> correct MSD) while pos_v reflects back in. Periodic keeps bmir=+1.
    double bmir[3] = { 1.0, 1.0, 1.0 };

    float* tr = traj + (size_t)w * (T + 1) * 3;
    tr[0] = pr[0]; tr[1] = pr[1]; tr[2] = pr[2];
    size_t tb = (size_t)w * (T + 1); t2sum[tb] = 0.0f; float t2acc = 0.0f;  // Phase 6b

    for (unsigned t = 1; t <= T; ++t) {
        // per-timestep step length from the CURRENT compartment (permeable crossings may
        // have changed in_obs).
        double step_length = (in_obs >= 0) ? obs_step[in_obs] : step_extra;
        // EPS just skips the ~0 root of the surface we currently rest on. Triangles are
        // single-root, so a blanket tc>EPS never tunnels to a "far root" (unlike a sphere);
        // it simply ignores the face the walker is sitting on (and any face sharing that
        // exact point, e.g. across a shared edge). Double-clean, << any real triangle gap.
        double EPS = 1e-9 * voxext;
        double gx = curand_normal(&st), gy = curand_normal(&st), gz = curand_normal(&st);
        double inv = 1.0 / sqrt(gx*gx + gy*gy + gz*gz + 1e-300);
        double dir[3] = { gx*inv, gy*inv, gz*inv };
        double remaining = step_length;

        for (int b = 0; b < 24 && remaining > EPS; ++b) {
            double O[3] = { pv[0], pv[1], pv[2] };
            double E[3] = { pv[0]+remaining*dir[0], pv[1]+remaining*dir[1], pv[2]+remaining*dir[2] };
            // grid cell range overlapping the step segment's AABB (clamped to the grid)
            int cmn[3], cmx[3];
            for (int a = 0; a < 3; ++a) {
                double lo = fmin(O[a], E[a]) - gmin[a], hi = fmax(O[a], E[a]) - gmin[a];
                int clo = int(floor(lo*ics[a])), chi = int(floor(hi*ics[a]));
                cmn[a] = clo < 0 ? 0 : (clo > dim[a]-1 ? dim[a]-1 : clo);
                cmx[a] = chi < 0 ? 0 : (chi > dim[a]-1 ? dim[a]-1 : chi);
            }
            double best_t = remaining; int best_tri = -1, best_axis = -1;
            for (int cz = cmn[2]; cz <= cmx[2]; ++cz)
            for (int cy = cmn[1]; cy <= cmx[1]; ++cy)
            for (int cx = cmn[0]; cx <= cmx[0]; ++cx) {
                int cell = cx + dim[0]*(cy + dim[1]*cz);
                for (int k = cell_start[cell]; k < cell_start[cell+1]; ++k) {
                    int i = cell_tri[k];
                    int i0 = tris[3*i+0], i1 = tris[3*i+1], i2 = tris[3*i+2];
                    double a[3] = { verts[3*i0+0], verts[3*i0+1], verts[3*i0+2] };
                    double bb[3]= { verts[3*i1+0], verts[3*i1+1], verts[3*i1+2] };
                    double cc[3]= { verts[3*i2+0], verts[3*i2+1], verts[3*i2+2] };
                    double tc;
                    if (mtHit(O, dir, a, bb, cc, tc) && tc > EPS && tc < best_t) { best_t = tc; best_tri = i; best_axis = -1; }
                }
            }
            for (int i = 0; i < 3; ++i) {               // voxel faces (in pv)
                double tw = -1.0;
                if (dir[i] >  1e-30) tw = (vmax[i]-pv[i])/dir[i];
                else if (dir[i] < -1e-30) tw = (vmin[i]-pv[i])/dir[i];
                if (tw > EPS && tw < best_t) { best_t = tw; best_axis = i; best_tri = -1; }
            }
            for (int i = 0; i < 3; ++i) { pr[i] += best_t*dir[i]*bmir[i]; pv[i] += best_t*dir[i]; }
            remaining -= best_t;
            if (best_tri >= 0) {                         // triangle membrane: cross or reflect
                int obs = tri_obs[best_tri];             // obstacle owning this surface
                bool from_intra = (in_obs == obs);       // inside obs -> i->e; else entering -> e->i
                float prob = from_intra ? prob_ie[obs] : prob_ei[obs];
                if (prob > 0.0f && prob > curand_uniform(&st)) {
                    // PERMEATE (CPU model): the membrane is transparent -- finish the step
                    // straight through (ignore further membranes) and re-classify at landing.
                    for (int i = 0; i < 3; ++i) { pr[i] += remaining*dir[i]*bmir[i]; pv[i] += remaining*dir[i]; }
                    remaining = 0.0;
                    in_obs = classifyMesh(pv, verts, tris, tri_obs, gmin, ics, dim, cell_start, cell_tri);
                } else {                                  // reflect off the triangle (double => lands exactly on the surface; no nudge)
                    int i0 = tris[3*best_tri+0], i1 = tris[3*best_tri+1], i2 = tris[3*best_tri+2];
                    double e1[3] = { verts[3*i1+0]-verts[3*i0+0], verts[3*i1+1]-verts[3*i0+1], verts[3*i1+2]-verts[3*i0+2] };
                    double e2[3] = { verts[3*i2+0]-verts[3*i0+0], verts[3*i2+1]-verts[3*i0+1], verts[3*i2+2]-verts[3*i0+2] };
                    double n[3] = { e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0] };
                    double nn = 1.0 / sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] + 1e-300);
                    n[0]*=nn; n[1]*=nn; n[2]*=nn;
                    double dn = dir[0]*n[0] + dir[1]*n[1] + dir[2]*n[2];
                    dir[0] -= 2.0*dn*n[0]; dir[1] -= 2.0*dn*n[1]; dir[2] -= 2.0*dn*n[2];
                }
            } else if (best_axis >= 0) {                 // voxel face
                if (periodic) {
                    pv[best_axis] = (dir[best_axis] > 0.0) ? vmin[best_axis] : vmax[best_axis];
                    // pr is NOT wrapped (continuous); bmir stays +1 -> pr already advanced straight.
                } else {
                    // MIRROR: flip dir AND bmir on this axis. pv reflects back into the voxel;
                    // pr's per-axis product dir*bmir is unchanged -> it keeps going straight (free).
                    dir[best_axis]  = -dir[best_axis];
                    bmir[best_axis] = -bmir[best_axis];
                }
            } else {
                remaining = 0.0;
            }
        }
        float T2c = (in_obs >= 0) ? obs_T2[in_obs] : t2_extra;  // T2 of the END-of-step compartment
        t2acc += dt / T2c; t2sum[tb + t] = t2acc;               // Phase 3b/6b: per-step T2 accumulation
        float* o = tr + (size_t)t * 3;
        o[0] = pr[0]; o[1] = pr[1]; o[2] = pr[2];
    }
    // Final illegal-crossing check (periodic/mirror-aware; mirrors CPU finalPositionCheck).
    // Re-classify the landing pos_v and compare its intra/extra MEMBERSHIP to the TRACKED
    // in_obs (which already follows every legitimate permeation). For a tiled periodic/mirror
    // substrate a wrap/fold preserves membership, so this fires ONLY on a real leak (the walker
    // tunnelled an impermeable membrane without in_obs being updated) -> flag -2 (signalKernel
    // drops it). Comparing MEMBERSHIP (not the obstacle index) is what makes it wrap/fold-safe.
    int final_obs = classifyMesh(pv, verts, tris, tri_obs, gmin, ics, dim, cell_start, cell_tri);
    out_inobs[w] = ((final_obs >= 0) != (in_obs >= 0)) ? -2 : final_obs;
}

//! Wrap a coordinate into [lo,hi) for PERIODIC BC: pv = lo + mod(x-lo, L) (sawtooth). For
//! periodic voxel BC this reproduces MCDC's incremental pv exactly (each face crossing sets
//! pv to the opposite face).
__device__ inline float wrapCoord(float x, float lo, float hi)
{
    float L = hi - lo;
    if (L <= 0.0f) return x;
    float t = x - lo;
    t = t - L * floorf(t / L);
    return lo + t;
}

//! Fold a coordinate into [lo,hi] for MIRROR BC: triangle wave (period 2L, reflect in the
//! upper half). pos_r runs free, so the in-voxel image pos_v = fold(pos_r) -- this reproduces
//! MCDC's incremental boundary_mirror pv exactly (each face crossing reflects pv back in).
__device__ inline float mirrorFold(float x, float lo, float hi)
{
    float L = hi - lo;
    if (L <= 0.0f) return x;
    float t = x - lo;
    float p = t - 2.0f*L * floorf(t / (2.0f*L));   // p in [0, 2L)
    return lo + (p <= L ? p : 2.0f*L - p);
}

//! Fold a real position into the voxel (sawtooth periodic / triangle-wave mirror) and return
//! the containing subdivision box index (-1 = outside). Regular grid (sub_n>0) -> O(1) index;
//! else an O(n_sub) scan over subdiv. Used for both final-position and per-echo-time binning.
__device__ inline int subBoxIndex(float px, float py, float pz,
                                  float vminx, float vminy, float vminz,
                                  float vmaxx, float vmaxy, float vmaxz, int periodic,
                                  const float* subdiv, unsigned n_sub, int sub_n)
{
    float pv[3] = { px, py, pz };
    if (periodic) {
        pv[0] = wrapCoord(pv[0], vminx, vmaxx); pv[1] = wrapCoord(pv[1], vminy, vmaxy); pv[2] = wrapCoord(pv[2], vminz, vmaxz);
    } else {
        pv[0] = mirrorFold(pv[0], vminx, vmaxx); pv[1] = mirrorFold(pv[1], vminy, vmaxy); pv[2] = mirrorFold(pv[2], vminz, vmaxz);
    }
    if (sub_n > 0) {
        float lo[3] = { vminx, vminy, vminz }, hi[3] = { vmaxx, vmaxy, vmaxz };
        int ijk[3]; bool inside = true;
        for (int a = 0; a < 3; ++a) {
            float gap = (hi[a] - lo[a]) / float(sub_n);
            int i = int(floorf((pv[a] - lo[a]) / gap));
            if (pv[a] < lo[a] || pv[a] > hi[a]) inside = false;
            if (i < 0) i = 0; else if (i > sub_n - 1) i = sub_n - 1;
            ijk[a] = i;
        }
        return inside ? (ijk[0] * sub_n + ijk[1]) * sub_n + ijk[2] : -1;
    }
    for (unsigned k = 0; k < n_sub; ++k) {
        const float* b = subdiv + (size_t)k * 6;
        if (pv[0] >= b[0] && pv[1] >= b[1] && pv[2] >= b[2] &&
            pv[0] <= b[3] && pv[1] <= b[4] && pv[2] <= b[5]) return int(k);
    }
    return -1;
}

//! One block per walker, num_rep measurements strided across threads. Reads the
//! on-device trajectory and accumulates the (T2=inf) signal sum_w cos(phase).
//! Mirrors CUDASynth::deviceDataSynth but uses the scheme values directly (MCDC).
// out_DWI layout: [0, 2*num_rep) is the GLOBAL signal -- [0,num_rep) intra (in_obs>=0),
// [num_rep,2*num_rep) extra. When n_sub>0, [2*num_rep ...) holds the per-subdivision
// VOLUME: subdivision k occupies [ (1+k)*2*num_rep, (2+k)*2*num_rep ) with the same
// intra|extra halves. A walker is binned by its FINAL wrapped position into the first
// subdivision box that contains it (mutually exclusive, matching the CPU's break-on-first),
// exactly as PGSESequence::update_DWI_signal does with walker.pos_v at echo time.
__global__ void signalKernel(const float* traj, unsigned T, float duration,
                             const float* scheme, unsigned num_rep,
                             const int* inobs, const float* t2sum,
                             float vminx, float vminy, float vminz,
                             float vmaxx, float vmaxy, float vmaxz, int periodic,
                             const float* subdiv, unsigned n_sub, int sub_n, int at_te, double* out_DWI)
{
    unsigned w = blockIdx.x;
    if (inobs[w] == -2) return;   // walker flagged as an illegal crossing (numerical membrane
                                  // leak) by the kernel's final re-classification -> excluded
                                  // from the signal, mirroring the CPU deportation. One block
                                  // per walker, so all threads return together (no __syncthreads).
    float time_step = duration / float(T);
    float dos_pi = 2.0f * float(M_PI);
    unsigned base = (inobs[w] >= 0) ? 0u : num_rep;   // intra -> first half, extra -> second
    const float* t2w_acc = t2sum + (size_t)w * (T + 1);  // per-step cumulative dt/T2 (Phase 6b)

    const float* tr = traj + (size_t)w * (T + 1) * 3;
    float x_ini[3] = { tr[0], tr[1], tr[2] };

    // Subdivision membership. Default (at_te==0): bin every acquisition by the FINAL position --
    // computed once per block in thread 0 and shared. With at_te==1 each acquisition is binned by
    // the walker's position AT ITS echo time (computed per-acquisition in the loop) -> density-vs-TE.
    __shared__ int sub_idx;
    if (n_sub > 0 && !at_te && threadIdx.x == 0) {
        const float* fp = tr + (size_t)T * 3;       // final pos_r (free / unwrapped)
        sub_idx = subBoxIndex(fp[0], fp[1], fp[2], vminx, vminy, vminz,
                              vmaxx, vmaxy, vmaxz, periodic, subdiv, n_sub, sub_n);
    }
    if (n_sub > 0 && !at_te) __syncthreads();

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
        // T2 relaxation weight at this measurement's echo time. t2sum[pos] is the
        // per-step accumulated dt/T2 over the walker's ACTUAL compartment path (handles
        // permeable crossings); = pos*dt/T2 for a fixed compartment. pos=min(int(TE/dt),T-1).
        int pos = int(te / time_step); if (pos >= int(T)) pos = int(T) - 1;
        float t2w = __expf(-t2w_acc[pos]);
        // Accumulate in DOUBLE: the global bin sums ~1 per walker over up to 1e8+ walkers
        // (~1e8), which a float32 accumulator cannot hold (24-bit mantissa saturates near
        // 2^24 ~ 1.7e7, so adds of ~1 are lost). double atomicAdd keeps it exact at scale.
        double contrib = double(t2w) * double(cosf(phase));
        atomicAdd(&out_DWI[base + s], contrib);                       // global signal
        if (n_sub > 0) {                                              // volumetric: this walker's box
            int sidx = sub_idx;                                       // final-position box (at_te==0)
            if (at_te) {                                             // per-echo-time: box at TE(s)
                const float* fp = tr + (size_t)pos * 3;             // pos_r at this acquisition's TE
                sidx = subBoxIndex(fp[0], fp[1], fp[2], vminx, vminy, vminz,
                                   vmaxx, vmaxy, vmaxz, periodic, subdiv, n_sub, sub_n);
            }
            if (sidx >= 0)
                atomicAdd(&out_DWI[(size_t)(1 + sidx) * 2u * num_rep + base + s], contrib);
        }
    }
}

static int gpuSynthAndFree(const GpuSimParams& p, float* d_traj, const int* d_inobs,
                           const float* d_t2sum, const float* scheme, double* out_DWI);   // defined below

int runGpuFreeDiffusion(const GpuSimParams& p, const float* scheme, double* out_DWI)
{
    float *d_traj = nullptr, *d_t2sum = nullptr; int *d_inobs = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj,  (size_t)p.N * (p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_inobs, (size_t)p.N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_t2sum, (size_t)p.N * (p.T + 1) * sizeof(float)));
    const float dt = p.duration / float(p.T);

    const int bs = 256;
    const int gs = int((p.N + bs - 1) / bs);
    freeDiffusionKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_length,
                                    p.vmin[0], p.vmin[1], p.vmin[2],
                                    p.vmax[0], p.vmax[1], p.vmax[2],
                                    p.smin[0], p.smin[1], p.smin[2],
                                    p.smax[0], p.smax[1], p.smax[2],
                                    p.periodic, dt, p.t2_extra, d_inobs, d_t2sum, p.seed);
    CUDA_CHECK(cudaGetLastError());

    int rc = gpuSynthAndFree(p, d_traj, d_inobs, d_t2sum, scheme, out_DWI);
    cudaFree(d_traj);
    cudaFree(d_inobs);
    cudaFree(d_t2sum);
    return rc;
}

// Shared tail: upload scheme, run the signal kernel on the (already-filled) device
// trajectory + per-walker compartment array, copy back the 2*num_rep DWI (intra half
// then extra half), free. out_DWI must hold 2*num_rep entries.
static int gpuSynthAndFree(const GpuSimParams& p, float* d_traj, const int* d_inobs,
                           const float* d_t2sum, const float* scheme, double* out_DWI)
{
    const size_t scheme_bytes = (size_t)p.num_rep * 7 * sizeof(float);
    // DWI buffer: global (2*num_rep) + one 2*num_rep block per subdivision (volumetric).
    const size_t n_blocks  = (size_t)1 + p.n_sub;
    const size_t dwi_count = n_blocks * p.num_rep * 2;
    const size_t dwi_bytes = dwi_count * sizeof(double);          // double accumulator (see signalKernel)
    float *d_scheme = nullptr, *d_subdiv = nullptr; double *d_dwi = nullptr;
    CUDA_CHECK(cudaMalloc(&d_scheme, scheme_bytes));
    CUDA_CHECK(cudaMalloc(&d_dwi,    dwi_bytes));
    CUDA_CHECK(cudaMemcpy(d_scheme, scheme, scheme_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_dwi, 0, dwi_bytes));
    if (p.n_sub > 0 && p.sub_n <= 0) {   // arbitrary boxes (file): upload them for the scan
        CUDA_CHECK(cudaMalloc(&d_subdiv, (size_t)p.n_sub * 6 * sizeof(float)));
        CUDA_CHECK(cudaMemcpy(d_subdiv, p.subdiv_boxes, (size_t)p.n_sub * 6 * sizeof(float),
                              cudaMemcpyHostToDevice));
    }

    const int threads = p.num_rep < 1024u ? int(p.num_rep) : 1024;
    signalKernel<<<p.N, threads>>>(d_traj, p.T, p.duration, d_scheme, p.num_rep,
                                   d_inobs, d_t2sum,
                                   p.vmin[0], p.vmin[1], p.vmin[2],
                                   p.vmax[0], p.vmax[1], p.vmax[2], p.periodic,
                                   d_subdiv, p.n_sub, p.sub_n, p.at_te, d_dwi);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(out_DWI, d_dwi, dwi_bytes, cudaMemcpyDeviceToHost));   // double -> double

    cudaFree(d_scheme);
    cudaFree(d_dwi);
    if (d_subdiv) cudaFree(d_subdiv);
    return 0;
}

int runGpuSpheres(const GpuSimParams& p, const float* spheres, const float* obs_T2,
                  const float* prob_ie, const float* prob_ei,
                  unsigned n_sph, const float* scheme, double* out_DWI)
{
    // Walker-chunked (same scheme as runGpuMesh): resident sphere data + persistent DWI
    // accumulator built once; only the per-chunk trajectory scales with N. Lets analytic
    // sphere runs reach large N (e.g. with a 1e6-subdivision accumulator) on a fixed GPU.
    float *d_sph = nullptr, *d_T2 = nullptr, *d_pie = nullptr, *d_pei = nullptr;
    CUDA_CHECK(cudaMalloc(&d_sph,  (size_t)n_sph * 4 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_T2,   (size_t)n_sph * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pie,  (size_t)n_sph * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pei,  (size_t)n_sph * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_sph, spheres, (size_t)n_sph * 4 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_T2,  obs_T2,  (size_t)n_sph * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pie, prob_ie, (size_t)n_sph * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pei, prob_ei, (size_t)n_sph * sizeof(float), cudaMemcpyHostToDevice));

    // --- uniform spatial grid over the voxel: CSR cell -> sphere-index lists, so each bounce
    // tests only spheres in the cells its step segment overlaps (O(nearby) vs O(n_sph)). Each
    // sphere is inserted into every cell its AABB (center +/- radius) touches; ~cbrt(n_sph) cells
    // per axis (=> ~1 cell for tiny lists, i.e. brute-force, and real pruning for dense gamma). ---
    double gmin[3], gmax[3];
    for (int a = 0; a < 3; ++a) { gmin[a] = p.vmin[a]; gmax[a] = p.vmax[a]; }
    int gper = int(cbrt(double(n_sph > 0 ? n_sph : 1))); if (gper < 1) gper = 1; if (gper > 128) gper = 128;
    int gdim[3]; double gics[3];
    for (int a = 0; a < 3; ++a) {
        double ext = gmax[a] - gmin[a]; double pad = 1e-4*ext + 1e-9;
        gmin[a] -= pad; gmax[a] += pad; ext = gmax[a] - gmin[a];
        gdim[a] = gper; gics[a] = double(gdim[a]) / ext;
    }
    const size_t g_ncells = size_t(gdim[0]) * gdim[1] * gdim[2];
    auto sphCellRange = [&](unsigned i, int cmn[3], int cmx[3]) {
        double r = spheres[i*4+3];
        for (int a = 0; a < 3; ++a) {
            double lo = (spheres[i*4+a] - r) - gmin[a], hi = (spheres[i*4+a] + r) - gmin[a];
            int clo = int(floor(lo*gics[a])), chi = int(floor(hi*gics[a]));
            cmn[a] = clo < 0 ? 0 : (clo > gdim[a]-1 ? gdim[a]-1 : clo);
            cmx[a] = chi < 0 ? 0 : (chi > gdim[a]-1 ? gdim[a]-1 : chi);
        }
    };
    std::vector<int> g_cstart(g_ncells + 1, 0); int cmn[3], cmx[3];
    for (unsigned i = 0; i < n_sph; ++i) { sphCellRange(i, cmn, cmx);
        for (int cz = cmn[2]; cz <= cmx[2]; ++cz) for (int cy = cmn[1]; cy <= cmx[1]; ++cy) for (int cx = cmn[0]; cx <= cmx[0]; ++cx)
            g_cstart[size_t(cx + gdim[0]*(cy + gdim[1]*cz)) + 1]++;
    }
    for (size_t c = 0; c < g_ncells; ++c) g_cstart[c+1] += g_cstart[c];
    std::vector<int> g_csph(g_cstart[g_ncells]);
    std::vector<int> g_cursor(g_cstart.begin(), g_cstart.end() - 1);
    for (unsigned i = 0; i < n_sph; ++i) { sphCellRange(i, cmn, cmx);
        for (int cz = cmn[2]; cz <= cmx[2]; ++cz) for (int cy = cmn[1]; cy <= cmx[1]; ++cy) for (int cx = cmn[0]; cx <= cmx[0]; ++cx)
            g_csph[g_cursor[size_t(cx + gdim[0]*(cy + gdim[1]*cz))]++] = int(i);
    }
    int *d_cstart = nullptr, *d_csph = nullptr;
    CUDA_CHECK(cudaMalloc(&d_cstart, (g_ncells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_csph,   (g_csph.empty() ? 1 : g_csph.size()) * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cstart, g_cstart.data(), (g_ncells + 1) * sizeof(int), cudaMemcpyHostToDevice));
    if (!g_csph.empty())
        CUDA_CHECK(cudaMemcpy(d_csph, g_csph.data(), g_csph.size() * sizeof(int), cudaMemcpyHostToDevice));

    const float dt = p.duration / float(p.T);

    const size_t dwi_count = ((size_t)1 + p.n_sub) * p.num_rep * 2;
    float *d_scheme = nullptr, *d_subdiv = nullptr; double *d_dwi = nullptr;
    CUDA_CHECK(cudaMalloc(&d_scheme, (size_t)p.num_rep * 7 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_dwi,    dwi_count * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_scheme, scheme, (size_t)p.num_rep * 7 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_dwi, 0, dwi_count * sizeof(double)));
    if (p.n_sub > 0 && p.sub_n <= 0) {
        CUDA_CHECK(cudaMalloc(&d_subdiv, (size_t)p.n_sub * 6 * sizeof(float)));
        CUDA_CHECK(cudaMemcpy(d_subdiv, p.subdiv_boxes, (size_t)p.n_sub * 6 * sizeof(float), cudaMemcpyHostToDevice));
    }

    const int bs = 256;
    const size_t per_walker = (size_t)(p.T + 1) * 3 * sizeof(float) + (size_t)(p.T + 1) * sizeof(float) + sizeof(int);
    size_t freeB = 0, totalB = 0;
    CUDA_CHECK(cudaMemGetInfo(&freeB, &totalB));
    size_t chunk = size_t(double(freeB) * 0.80) / per_walker;
    if (chunk >= (size_t)p.N) chunk = p.N;
    else { chunk = (chunk / bs) * bs; if (chunk < (size_t)bs) { std::fprintf(stderr, "[GPU][sph][ERROR] not enough VRAM.\n"); return int(cudaErrorMemoryAllocation); } }
    if (const char* e = std::getenv("MCDC_GPU_CHUNK")) { long long c = std::atoll(e); if (c > 0 && (size_t)c < chunk) chunk = (size_t)c; }
    const unsigned n_chunks = unsigned((p.N + chunk - 1) / chunk);
    std::fprintf(stderr, "[GPU][sph] N=%u, T=%u: free=%.2f/%.2f GB -> chunk=%zu walkers x %u chunk(s).\n",
                 p.N, p.T, freeB/1e9, totalB/1e9, chunk, n_chunks);

    float *d_traj = nullptr, *d_t2sum = nullptr; int *d_inobs = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj,  chunk * (size_t)(p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_t2sum, chunk * (size_t)(p.T + 1) * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_inobs, chunk * sizeof(int)));
    const int sig_threads = p.num_rep < 1024u ? int(p.num_rep) : 1024;

    int rc = 0;
    auto t_start = std::chrono::steady_clock::now(); int last_pct = -1;
    for (unsigned c0 = 0; c0 < p.N; c0 += unsigned(chunk)) {
        unsigned cn = (p.N - c0 < chunk) ? (p.N - c0) : unsigned(chunk);
        int gs = int((cn + bs - 1) / bs);
        spheresKernel<<<gs, bs>>>(d_traj, cn, p.T, p.step_intra, p.step_extra, d_sph, n_sph,
                                  d_pie, d_pei,
                                  d_cstart, d_csph,
                                  gmin[0], gmin[1], gmin[2], gics[0], gics[1], gics[2],
                                  gdim[0], gdim[1], gdim[2],
                                  p.vmin[0], p.vmin[1], p.vmin[2], p.vmax[0], p.vmax[1], p.vmax[2],
                                  p.smin[0], p.smin[1], p.smin[2], p.smax[0], p.smax[1], p.smax[2],
                                  p.periodic, p.placement, dt, p.t2_extra, d_T2, d_inobs, d_t2sum, p.seed, c0);
        if ((rc = int(cudaGetLastError()))) break;
        signalKernel<<<cn, sig_threads>>>(d_traj, p.T, p.duration, d_scheme, p.num_rep,
                                          d_inobs, d_t2sum,
                                          p.vmin[0], p.vmin[1], p.vmin[2], p.vmax[0], p.vmax[1], p.vmax[2], p.periodic,
                                          d_subdiv, p.n_sub, p.sub_n, p.at_te, d_dwi);
        if ((rc = int(cudaGetLastError()))) break;
        if ((rc = int(cudaDeviceSynchronize()))) break;
        unsigned done_w = c0 + cn; int pct = int(100.0 * double(done_w) / double(p.N));
        if (pct != last_pct || c0 == 0) {
            last_pct = pct;
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
            double rate = double(done_w) / (el > 0 ? el : 1e-9);
            std::fprintf(stderr, "[GPU][sph] %3d%% | %.0fs | ETA %.0fs | %.0f walkers/s\n",
                         pct, el, (rate > 0 ? double(p.N - done_w)/rate : 0.0), rate);
        }
    }
    if (rc == 0) CUDA_CHECK(cudaMemcpy(out_DWI, d_dwi, dwi_count * sizeof(double), cudaMemcpyDeviceToHost));

    cudaFree(d_traj); cudaFree(d_t2sum); cudaFree(d_inobs);
    cudaFree(d_scheme); cudaFree(d_dwi); if (d_subdiv) cudaFree(d_subdiv);
    cudaFree(d_sph); cudaFree(d_T2); cudaFree(d_pie); cudaFree(d_pei);
    cudaFree(d_cstart); cudaFree(d_csph);
    return rc;
}

int runGpuCylinders(const GpuSimParams& p, const float* cyls, const float* obs_T2,
                    const float* prob_ie, const float* prob_ei,
                    unsigned n_cyl, const float* scheme, double* out_DWI)
{
    float *d_traj = nullptr, *d_cyl = nullptr, *d_T2 = nullptr, *d_pie = nullptr, *d_pei = nullptr, *d_t2sum = nullptr;
    int *d_inobs = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj, (size_t)p.N * (p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_cyl,  (size_t)n_cyl * 7 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_T2,   (size_t)n_cyl * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pie,  (size_t)n_cyl * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pei,  (size_t)n_cyl * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_inobs,(size_t)p.N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_t2sum,(size_t)p.N * (p.T + 1) * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_cyl, cyls, (size_t)n_cyl * 7 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_T2,  obs_T2,  (size_t)n_cyl * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pie, prob_ie, (size_t)n_cyl * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pei, prob_ei, (size_t)n_cyl * sizeof(float), cudaMemcpyHostToDevice));

    // --- spatial grid over the voxel: CSR cell -> cylinder-index lists (O(nearby) collision).
    // A cylinder is FINITE along a world axis only where its direction component is ~0; along a
    // direction it is parallel to, it spans the whole grid. So axis-aligned (z-parallel) packings
    // prune in the perpendicular plane; tilted cylinders span more cells (=> ~brute, still correct). ---
    double gmin[3], gmax[3];
    for (int a = 0; a < 3; ++a) { gmin[a] = p.vmin[a]; gmax[a] = p.vmax[a]; }
    int gper = int(cbrt(double(n_cyl > 0 ? n_cyl : 1))); if (gper < 1) gper = 1; if (gper > 128) gper = 128;
    int gdim[3]; double gics[3];
    for (int a = 0; a < 3; ++a) {
        double ext = gmax[a] - gmin[a]; double pad = 1e-4*ext + 1e-9;
        gmin[a] -= pad; gmax[a] += pad; ext = gmax[a] - gmin[a];
        gdim[a] = gper; gics[a] = double(gdim[a]) / ext;
    }
    const size_t g_ncells = size_t(gdim[0]) * gdim[1] * gdim[2];
    auto cylCellRange = [&](unsigned i, int cmn[3], int cmx[3]) {
        double r = cyls[i*7+6];
        for (int a = 0; a < 3; ++a) {
            double Da = cyls[i*7+3+a];
            double loW, hiW;
            if (fabs(Da) < 1e-6) { loW = cyls[i*7+a] - r; hiW = cyls[i*7+a] + r; }  // finite along a
            else                 { loW = gmin[a];         hiW = gmax[a];         }  // parallel -> span grid
            int clo = int(floor((loW - gmin[a])*gics[a])), chi = int(floor((hiW - gmin[a])*gics[a]));
            cmn[a] = clo < 0 ? 0 : (clo > gdim[a]-1 ? gdim[a]-1 : clo);
            cmx[a] = chi < 0 ? 0 : (chi > gdim[a]-1 ? gdim[a]-1 : chi);
        }
    };
    std::vector<int> g_cstart(g_ncells + 1, 0); int cmn[3], cmx[3];
    for (unsigned i = 0; i < n_cyl; ++i) { cylCellRange(i, cmn, cmx);
        for (int cz = cmn[2]; cz <= cmx[2]; ++cz) for (int cy = cmn[1]; cy <= cmx[1]; ++cy) for (int cx = cmn[0]; cx <= cmx[0]; ++cx)
            g_cstart[size_t(cx + gdim[0]*(cy + gdim[1]*cz)) + 1]++;
    }
    for (size_t c = 0; c < g_ncells; ++c) g_cstart[c+1] += g_cstart[c];
    std::vector<int> g_ccyl(g_cstart[g_ncells]);
    std::vector<int> g_cursor(g_cstart.begin(), g_cstart.end() - 1);
    for (unsigned i = 0; i < n_cyl; ++i) { cylCellRange(i, cmn, cmx);
        for (int cz = cmn[2]; cz <= cmx[2]; ++cz) for (int cy = cmn[1]; cy <= cmx[1]; ++cy) for (int cx = cmn[0]; cx <= cmx[0]; ++cx)
            g_ccyl[g_cursor[size_t(cx + gdim[0]*(cy + gdim[1]*cz))]++] = int(i);
    }
    int *d_cstart = nullptr, *d_ccyl = nullptr;
    CUDA_CHECK(cudaMalloc(&d_cstart, (g_ncells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_ccyl,   (g_ccyl.empty() ? 1 : g_ccyl.size()) * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_cstart, g_cstart.data(), (g_ncells + 1) * sizeof(int), cudaMemcpyHostToDevice));
    if (!g_ccyl.empty())
        CUDA_CHECK(cudaMemcpy(d_ccyl, g_ccyl.data(), g_ccyl.size() * sizeof(int), cudaMemcpyHostToDevice));

    const float dt = p.duration / float(p.T);

    const int bs = 256, gs = int((p.N + bs - 1) / bs);
    cylindersKernel<<<gs, bs>>>(d_traj, p.N, p.T, p.step_intra, p.step_extra, d_cyl, n_cyl,
                                d_pie, d_pei,
                                d_cstart, d_ccyl,
                                gmin[0], gmin[1], gmin[2], gics[0], gics[1], gics[2],
                                gdim[0], gdim[1], gdim[2],
                                p.vmin[0], p.vmin[1], p.vmin[2], p.vmax[0], p.vmax[1], p.vmax[2],
                                p.smin[0], p.smin[1], p.smin[2], p.smax[0], p.smax[1], p.smax[2],
                                p.periodic, p.placement, dt, p.t2_extra, d_T2, d_inobs, d_t2sum, p.seed);
    CUDA_CHECK(cudaGetLastError());

    int rc = gpuSynthAndFree(p, d_traj, d_inobs, d_t2sum, scheme, out_DWI);
    cudaFree(d_traj); cudaFree(d_cyl); cudaFree(d_T2);
    cudaFree(d_pie);  cudaFree(d_pei); cudaFree(d_inobs); cudaFree(d_t2sum);
    cudaFree(d_cstart); cudaFree(d_ccyl);
    return rc;
}

// Host: cell range (inclusive, clamped) overlapping triangle i's vertex AABB.
static inline void triCellRange(const float* verts, const int* tris, unsigned i,
                                const float gmin[3], const float ics[3], const int dim[3],
                                int cmin[3], int cmax[3])
{
    int i0 = tris[3*i+0], i1 = tris[3*i+1], i2 = tris[3*i+2];
    for (int a = 0; a < 3; ++a) {
        float v0 = verts[3*i0+a], v1 = verts[3*i1+a], v2 = verts[3*i2+a];
        float lo = fminf(v0, fminf(v1, v2)) - gmin[a];
        float hi = fmaxf(v0, fmaxf(v1, v2)) - gmin[a];
        int clo = int(floorf(lo*ics[a])), chi = int(floorf(hi*ics[a]));
        cmin[a] = clo < 0 ? 0 : (clo > dim[a]-1 ? dim[a]-1 : clo);
        cmax[a] = chi < 0 ? 0 : (chi > dim[a]-1 ? dim[a]-1 : chi);
    }
}

int runGpuMesh(const GpuSimParams& p, const float* verts, unsigned n_verts,
               const int* tris, unsigned n_tris,
               const int* tri_obs, unsigned n_obs,
               const float* obs_T2, const float* obs_step,
               const float* prob_ie, const float* prob_ei,
               const float vmin[3], const float vmax[3],
               const float* scheme, double* out_DWI)
{
    // --- build a uniform AABB grid (CSR cell->triangle lists) over the mesh bbox ---
    float gmin[3] = { 1e30f, 1e30f, 1e30f }, gmax[3] = { -1e30f, -1e30f, -1e30f };
    for (unsigned v = 0; v < n_verts; ++v)
        for (int a = 0; a < 3; ++a) {
            float c = verts[3*v+a];
            if (c < gmin[a]) gmin[a] = c;
            if (c > gmax[a]) gmax[a] = c;
        }
    int per = int(cbrtf(float(n_tris)));               // ~cubic-root cells per axis
    if (per < 1) per = 1; if (per > 256) per = 256;
    int dim[3]; float ics[3];
    for (int a = 0; a < 3; ++a) {
        float ext = gmax[a] - gmin[a];
        float pad = 1e-4f*ext + 1e-9f;                 // keep max-face points strictly inside
        gmin[a] -= pad; gmax[a] += pad; ext = gmax[a] - gmin[a];
        dim[a] = per;
        ics[a] = float(dim[a]) / ext;
    }
    const size_t ncells = size_t(dim[0]) * dim[1] * dim[2];

    std::vector<int> cell_start(ncells + 1, 0);
    int cmin[3], cmax[3];
    for (unsigned i = 0; i < n_tris; ++i) {            // pass 1: counts
        triCellRange(verts, tris, i, gmin, ics, dim, cmin, cmax);
        for (int cz = cmin[2]; cz <= cmax[2]; ++cz)
        for (int cy = cmin[1]; cy <= cmax[1]; ++cy)
        for (int cx = cmin[0]; cx <= cmax[0]; ++cx)
            cell_start[size_t(cx + dim[0]*(cy + dim[1]*cz)) + 1]++;
    }
    for (size_t c = 0; c < ncells; ++c) cell_start[c+1] += cell_start[c];
    std::vector<int> cell_tri(cell_start[ncells]);
    std::vector<int> cursor(cell_start.begin(), cell_start.end() - 1);
    for (unsigned i = 0; i < n_tris; ++i) {            // pass 2: fill
        triCellRange(verts, tris, i, gmin, ics, dim, cmin, cmax);
        for (int cz = cmin[2]; cz <= cmax[2]; ++cz)
        for (int cy = cmin[1]; cy <= cmax[1]; ++cy)
        for (int cx = cmin[0]; cx <= cmax[0]; ++cx)
            cell_tri[cursor[size_t(cx + dim[0]*(cy + dim[1]*cz))]++] = int(i);
    }

    // ---- resident device data: mesh + grid + per-obstacle arrays (built ONCE, reused
    // across every walker chunk). The trajectory is what blows up with N, so it is NOT
    // allocated here -- only per-chunk below. ----
    float *d_verts = nullptr, *d_oT2 = nullptr, *d_ostep = nullptr, *d_pie = nullptr, *d_pei = nullptr;
    int   *d_tris = nullptr, *d_cstart = nullptr, *d_ctri = nullptr, *d_tobs = nullptr;
    CUDA_CHECK(cudaMalloc(&d_verts,  (size_t)n_verts * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_tris,   (size_t)n_tris  * 3 * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_tobs,   (size_t)n_tris      * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_oT2,    (size_t)n_obs       * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_ostep,  (size_t)n_obs       * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pie,    (size_t)n_obs       * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_pei,    (size_t)n_obs       * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_cstart, (ncells + 1) * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_ctri,   (cell_tri.empty() ? 1 : cell_tri.size()) * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_verts,  verts, (size_t)n_verts * 3 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tris,   tris,  (size_t)n_tris  * 3 * sizeof(int),   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_tobs,   tri_obs, (size_t)n_tris    * sizeof(int),   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_oT2,    obs_T2,  (size_t)n_obs     * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ostep,  obs_step,(size_t)n_obs     * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pie,    prob_ie, (size_t)n_obs     * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pei,    prob_ei, (size_t)n_obs     * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cstart, cell_start.data(), (ncells + 1) * sizeof(int), cudaMemcpyHostToDevice));
    if (!cell_tri.empty())
        CUDA_CHECK(cudaMemcpy(d_ctri, cell_tri.data(), cell_tri.size() * sizeof(int), cudaMemcpyHostToDevice));
    const float dt = p.duration / float(p.T);

    // ---- persistent signal accumulator (spans all chunks): scheme + subdivisions + DWI.
    // DWI layout [global | per-subdivision blocks]; signalKernel atomicAdds into it. ----
    const size_t dwi_count = ((size_t)1 + p.n_sub) * p.num_rep * 2;
    float *d_scheme = nullptr, *d_subdiv = nullptr; double *d_dwi = nullptr;   // double accumulator
    CUDA_CHECK(cudaMalloc(&d_scheme, (size_t)p.num_rep * 7 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_dwi,    dwi_count * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_scheme, scheme, (size_t)p.num_rep * 7 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_dwi, 0, dwi_count * sizeof(double)));
    if (p.n_sub > 0 && p.sub_n <= 0) {   // arbitrary boxes (file): upload for the scan. A regular
        CUDA_CHECK(cudaMalloc(&d_subdiv, (size_t)p.n_sub * 6 * sizeof(float)));  // grid (sub_n>0) needs
        CUDA_CHECK(cudaMemcpy(d_subdiv, p.subdiv_boxes, (size_t)p.n_sub * 6 * sizeof(float), cudaMemcpyHostToDevice));
    }                                    // no boxes (O(1) arithmetic) -> saves n_sub*6 floats of VRAM

    // ---- size the walker CHUNK from the actually-free VRAM (after the resident data),
    // not "as much as possible": a safe fraction of free memory, rounded DOWN to whole
    // blocks so every launch is full thread-blocks (good occupancy, no ragged grids). ----
    const int bs = 256;
    const size_t per_walker = (size_t)(p.T + 1) * 3 * sizeof(float)   // trajectory (dominant)
                            + (size_t)(p.T + 1) * sizeof(float)       // t2sum
                            + sizeof(int);                            // in_obs
    size_t freeB = 0, totalB = 0;
    CUDA_CHECK(cudaMemGetInfo(&freeB, &totalB));
    size_t chunk = size_t(double(freeB) * 0.80) / per_walker;         // 20% VRAM headroom
    if (chunk >= (size_t)p.N) {
        chunk = p.N;                                                  // all fits -> single chunk (no rounding)
    } else {
        chunk = (chunk / bs) * bs;                                    // whole blocks per chunk
        if (chunk < (size_t)bs) {
            std::fprintf(stderr, "[GPU][mesh][ERROR] not enough VRAM: free=%.2f GB, need >=%zu B/walker.\n",
                         freeB/1e9, per_walker);
            return int(cudaErrorMemoryAllocation);
        }
    }
    // Optional override (testing / debugging): MCDC_GPU_CHUNK forces a smaller chunk so the
    // multi-chunk path can be exercised on small meshes. Results must be identical to the
    // single-chunk run (seeding is by global walker id), only the memory footprint changes.
    if (const char* e = std::getenv("MCDC_GPU_CHUNK")) {
        long long c = std::atoll(e);
        if (c > 0 && (size_t)c < chunk) chunk = (size_t)c;
    }
    const unsigned n_chunks = unsigned((p.N + chunk - 1) / chunk);
    std::fprintf(stderr, "[GPU][mesh] N=%u, T=%u: free=%.2f/%.2f GB, %zu B/walker -> chunk=%zu walkers x %u chunk(s), block=%d.\n",
                 p.N, p.T, freeB/1e9, totalB/1e9, per_walker, chunk, n_chunks, bs);

    // ---- per-chunk buffers allocated ONCE at the chunk size, reused every iteration ----
    float *d_traj = nullptr, *d_t2sum = nullptr; int *d_inobs = nullptr;
    CUDA_CHECK(cudaMalloc(&d_traj,  chunk * (size_t)(p.T + 1) * 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_t2sum, chunk * (size_t)(p.T + 1) * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_inobs, chunk * sizeof(int)));
    const int sig_threads = p.num_rep < 1024u ? int(p.num_rep) : 1024;

    int rc = 0;
    auto t_start = std::chrono::steady_clock::now();
    int last_pct = -1;
    for (unsigned c0 = 0; c0 < p.N; c0 += unsigned(chunk)) {
        unsigned ci = unsigned(c0 / chunk);                                // 0-based chunk index
        unsigned cn = (p.N - c0 < chunk) ? (p.N - c0) : unsigned(chunk);   // last chunk may be smaller
        int gs = int((cn + bs - 1) / bs);
        meshKernel<<<gs, bs>>>(d_traj, cn, p.T, p.step_extra, d_verts, d_tris, n_tris,
                               d_tobs, d_oT2, d_ostep, d_pie, d_pei,
                               vmin[0], vmin[1], vmin[2], vmax[0], vmax[1], vmax[2], p.periodic, p.placement,
                               p.smin[0], p.smin[1], p.smin[2], p.smax[0], p.smax[1], p.smax[2],
                               gmin[0], gmin[1], gmin[2], ics[0], ics[1], ics[2],
                               dim[0], dim[1], dim[2], d_cstart, d_ctri,
                               dt, p.t2_extra, d_inobs, d_t2sum, p.seed, c0);
        if ((rc = int(cudaGetLastError()))) break;
        signalKernel<<<cn, sig_threads>>>(d_traj, p.T, p.duration, d_scheme, p.num_rep,
                                          d_inobs, d_t2sum,
                                          vmin[0], vmin[1], vmin[2], vmax[0], vmax[1], vmax[2], p.periodic,
                                          d_subdiv, p.n_sub, p.sub_n, p.at_te, d_dwi);
        if ((rc = int(cudaGetLastError()))) break;
        if ((rc = int(cudaDeviceSynchronize()))) break;   // finish chunk before its buffers are reused

        // GPU progress + ETA, extrapolated from completed chunks (like the CPU progress bar).
        // Printed when the integer percentage advances, plus the very first chunk (immediate ETA).
        unsigned done_w = (c0 + cn);
        int pct = int(100.0 * double(done_w) / double(p.N));
        if (pct != last_pct || ci == 0) {
            last_pct = pct;
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
            double rate = double(done_w) / (el > 0 ? el : 1e-9);          // walkers/s so far
            double eta = (rate > 0) ? double(p.N - done_w) / rate : 0.0;  // seconds remaining
            std::fprintf(stderr, "[GPU][mesh] %3d%% | chunk %u/%u | %.0fs elapsed | ETA %.0fs (%.0fm) | %.0f walkers/s\n",
                         pct, ci + 1, n_chunks, el, eta, eta/60.0, rate);
        }
    }

    if (rc == 0) {
        CUDA_CHECK(cudaMemcpy(out_DWI, d_dwi, dwi_count * sizeof(double), cudaMemcpyDeviceToHost));
    }
    cudaFree(d_traj); cudaFree(d_t2sum); cudaFree(d_inobs);
    cudaFree(d_scheme); cudaFree(d_dwi); if (d_subdiv) cudaFree(d_subdiv);
    cudaFree(d_verts); cudaFree(d_tris); cudaFree(d_tobs);
    cudaFree(d_oT2); cudaFree(d_ostep); cudaFree(d_pie); cudaFree(d_pei);
    cudaFree(d_cstart); cudaFree(d_ctri);
    return rc;
}
