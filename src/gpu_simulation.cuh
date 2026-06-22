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
};

//! A single analytic sphere obstacle (internal units, mm).
struct GpuSphere
{
    float cx, cy, cz;   //!< centre
    float r;            //!< radius
};

//! \brief Free-diffusion dynamics + PGSE signal, fully on the GPU.
//! \param p        simulation parameters (internal units).
//! \param scheme   flattened scheme, num_rep*7 floats per row
//!                 (gx,gy,gz,G,Delta,delta,TE) exactly as stored by PGSESequence.
//! \param out_DWI  host array (size num_rep); filled with sum_w cos(phase) per row.
//! \return 0 on success, non-zero CUDA error code otherwise.
int runGpuFreeDiffusion(const GpuSimParams& p, const float* scheme, double* out_DWI);

//! \brief Phase 2a: restricted diffusion inside ONE impermeable sphere + signal.
//! Walkers start uniformly inside the sphere and elastically reflect off the inner
//! surface (brute-force, no acceleration grid). Same signal synthesis as above.
int runGpuSphere(const GpuSimParams& p, const GpuSphere& s, const float* scheme, double* out_DWI);

//! A single analytic (infinite) cylinder obstacle (internal units, mm).
struct GpuCylinder
{
    float px, py, pz;   //!< a point on the axis
    float dx, dy, dz;   //!< unit axis direction
    float r;            //!< radius
};

//! \brief Phase 2b: restricted diffusion inside ONE impermeable infinite cylinder.
//! Radially confined, free along the axis; reflects off the radial normal.
int runGpuCylinder(const GpuSimParams& p, const GpuCylinder& c, const float* scheme, double* out_DWI);

#endif // GPU_SIMULATION_CUH
