//!  GPU simulation interface (Phase 1: free diffusion + PGSE signal) =========
/*! \file gpu_simulation.cuh
 *  \brief Plain (Eigen-free) interface between the host C++ simulator and the CUDA
 *  device code. The host fills these POD structures from Parameters / PGSESequence
 *  (the "interface function"), so the .cu translation unit never sees Eigen or the
 *  heavy host classes. Phase 1 does free diffusion (no obstacles) and computes the
 *  PGSE signal on-device, reusing MCDC's exact internal-unit math (mm, ms, T).
 *  T2 is assumed infinite (no relaxation weighting) in this phase.
 *==============================================================================*/
#ifndef GPU_SIMULATION_CUH
#define GPU_SIMULATION_CUH

//! POD parameters for one GPU dynamics+synthesis launch. All in MCDC internal
//! units (lengths mm, time ms) to match the CPU path exactly.
struct GpuSimParams
{
    unsigned       N;            //!< number of walkers
    unsigned       T;            //!< number of timesteps
    unsigned       num_rep;      //!< number of measurements (scheme rows)
    float          duration;     //!< total dynamics duration (internal ms)
    float          step_length;  //!< sqrt(6*D*duration/T) (internal mm)
    unsigned long long seed;     //!< RNG seed
    float          vmin[3];      //!< voxel min corner (internal mm)
    float          vmax[3];      //!< voxel max corner (internal mm)
    int            periodic;     //!< voxel boundary: 1 = periodic, 0 = mirror
    int            placement;    //!< walker start: 0 = uniform in voxel, 1 = intra, 2 = extra
};

//! \brief Free-diffusion dynamics + PGSE signal, fully on the GPU.
//! \param p        simulation parameters (internal units).
//! \param scheme   flattened scheme, num_rep*7 floats per row
//!                 (gx,gy,gz,G,Delta,delta,TE) exactly as stored by PGSESequence.
//! \param out_DWI  host array (size num_rep); filled with sum_w cos(phase) per row.
//! \return 0 on success, non-zero CUDA error code otherwise.
int runGpuFreeDiffusion(const GpuSimParams& p, const float* scheme, double* out_DWI);

//! \brief Phase 2a/2c: restricted diffusion inside impermeable spheres + signal.
//! Walkers start uniformly inside a (random) sphere and reflect off the nearest
//! surface each step (brute-force over all spheres; no acceleration grid).
//! \param spheres  n_sph * 4 floats: (cx, cy, cz, r) per sphere (internal mm).
int runGpuSpheres(const GpuSimParams& p, const float* spheres, unsigned n_sph,
                  const float* scheme, double* out_DWI);

//! \brief Phase 2b/2c: restricted diffusion inside impermeable infinite cylinders.
//! Radially confined, free along each axis; reflects off the nearest radial normal.
//! \param cyls     n_cyl * 7 floats: (px,py,pz, dx,dy,dz, r) per cylinder (mm; D unit).
int runGpuCylinders(const GpuSimParams& p, const float* cyls, unsigned n_cyl,
                    const float* scheme, double* out_DWI);

//! \brief Phase 2d: restricted diffusion inside ONE triangle mesh (brute-force
//! Moller-Trumbore over all triangles; no acceleration grid). Walkers start
//! uniformly in [vmin,vmax] (assumed inside the closed mesh) and reflect off the
//! nearest triangle. \param verts 3*n_verts floats; \param tris 3*n_tris vertex indices.
int runGpuMesh(const GpuSimParams& p, const float* verts, unsigned n_verts,
               const int* tris, unsigned n_tris,
               const float vmin[3], const float vmax[3],
               const float* scheme, double* out_DWI);

#endif // GPU_SIMULATION_CUH
