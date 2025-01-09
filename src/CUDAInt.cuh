// CUDAInt.cuh
#ifndef CUDAINT_H_
#define CUDAINT_H_

#include "triangle.h"
#ifdef __CUDACC__
    // Device code only

    #include <cuda_runtime.h>  // For double3

    struct TriangleGPU
    {
        double3 A, B, C; // Each vertex: x, y, z
    };
    
    // Kernel Declaration
    __global__ void checkAllPairsKernel(const TriangleGPU* trianglesA, int sizeA,
                                    const TriangleGPU* trianglesB, int sizeB,
                                    int tileA_start, int tileB_start,
                                    int* result);
    
    // Device Function Declarations
    __device__ bool triangleIntersectsGPU(const double3 &A1,
                                          const double3 &B1,
                                          const double3 &C1,
                                          const double3 &A2,
                                          const double3 &B2,
                                          const double3 &C2);

    __device__ bool anyEdgeIntersectsTriangleGPU(const double3 &A1,
                                                 const double3 &B1,
                                                 const double3 &C1,
                                                 const double3 &A2,
                                                 const double3 &B2,
                                                 const double3 &C2);

    __device__ bool segmentIntersectsTriangleGPU(const double3 &P0,
                                                 const double3 &P1,
                                                 const double3 &A,
                                                 const double3 &B,
                                                 const double3 &C);

    __device__ bool isCloseGPU(const double3 &p1, const double3 &p2, double eps=1e-12);

    __global__ void checkGridBasedKernel(const TriangleGPU* triangles, int numTriangles,
                                     const int* gridCellOffsets, const int* gridCellTriangles,
                                     const int* gridDimensions, double3 minBounds,
                                     double3 cellSize, int* result);

    __global__ void checkGridBasedKernel(const TriangleGPU* triangles, int numTriangles,
                                        const int* gridCellOffsets, const int* gridCellTriangles,
                                        int3 gridDimensions, double3 minBounds,
                                        double3 cellSize, int* result);

    // Include device helper functions
    #include "device_utils.cuh"

#else
    // Host code only

    #include <string>
    #include <fstream>
    #include <vector>
    #include <iostream>
    #include "CUDAInt.cuh"  // To include any shared declarations if necessary

    // Host-side structures and function declarations can be added here
    // Example: Host function to run the collision kernel
    void runCollisionKernel(const std::vector<TriangleGPU>& h_triangles, std::vector<int>& h_result);

#endif // __CUDACC__

#endif /* CUDAINT_H_ */
