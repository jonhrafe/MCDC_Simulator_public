// CUDAInt.cu

#include "CUDAInt.cuh"
#include <cuda_runtime.h>
#include <cmath>

// Device Function Definitions

__device__ bool isCloseGPU(const double3 &p1, const double3 &p2, double eps) {
    // squared distance < eps^2
    double3 d = p1 - p2;
    double dist2 = dot3(d, d);
    return (dist2 < eps * eps);
}

// Segment-triangle intersection (Möller-Trumbore).
__device__ bool segmentIntersectsTriangleGPU(const double3 &P0,
                                             const double3 &P1,
                                             const double3 &A,
                                             const double3 &B,
                                             const double3 &C)
{
    const double EPS = 1e-10;
    double3 e1 = B - A;
    double3 e2 = C - A;
    double3 d  = P1 - P0;  // segment direction

    double3 p  = cross3(d, e2);
    double det = dot3(e1, p);

    if (fabs(det) < EPS) {
        return false; // parallel or no intersection
    }

    double invDet = 1.0 / det;
    double3 T = P0 - A;

    // Calculate U parameter
    double u = dot3(T, p) * invDet;
    if (u < 0.0 - EPS || u > 1.0 + EPS) {
        return false;
    }

    // Calculate V parameter
    double3 q = cross3(T, e1);
    double v = dot3(d, q) * invDet;
    if (v < 0.0 - EPS || (u + v) > 1.0 + EPS) {
        return false;
    }

    // Intersection "distance" along [P0,P1]
    double t = dot3(e2, q) * invDet;
    // For a segment, t must be in [0..1]
    if (t < 0.0 - EPS || t > 1.0 + EPS) {
        return false;
    }

    return true;
}

// Check if any edge of one triangle intersects the other triangle
__device__ bool anyEdgeIntersectsTriangleGPU(const double3 &A1,
                                             const double3 &B1,
                                             const double3 &C1,
                                             const double3 &A2,
                                             const double3 &B2,
                                             const double3 &C2)
{
    if (segmentIntersectsTriangleGPU(A1, B1, A2, B2, C2)) return true;
    if (segmentIntersectsTriangleGPU(B1, C1, A2, B2, C2)) return true;
    if (segmentIntersectsTriangleGPU(C1, A1, A2, B2, C2)) return true;
    return false;
}

// The main device function: does triangle-triangle intersection
__device__ bool triangleIntersectsGPU(const double3 &A1,
                                      const double3 &B1,
                                      const double3 &C1,
                                      const double3 &A2,
                                      const double3 &B2,
                                      const double3 &C2)
{
    // 1) Check adjacency (if they share a vertex, skip)
    for (int i = 0; i < 3; ++i) {
        double3 v1 = (i==0 ? A1 : (i==1 ? B1 : C1));
        for (int j = 0; j < 3; ++j) {
            double3 v2 = (j==0 ? A2 : (j==1 ? B2 : C2));
            if (isCloseGPU(v1, v2, 1e-12)) {
                return false; // adjacency => skip
            }
        }
    }

    // 2) Check edges
    if (anyEdgeIntersectsTriangleGPU(A1, B1, C1, A2, B2, C2)) return true;
    if (anyEdgeIntersectsTriangleGPU(A2, B2, C2, A1, B1, C1)) return true;

    return false;
}

// Kernel: check all pairs (i, j) of triangles
__global__ void checkAllPairsKernel(const TriangleGPU* trianglesA, int sizeA,
                                    const TriangleGPU* trianglesB, int sizeB,
                                    int tileA_start, int tileB_start,
                                    int* result)
{
    int local_i = blockIdx.x * blockDim.x + threadIdx.x;
    int local_j = blockIdx.y * blockDim.y + threadIdx.y;

    if (local_i < sizeA && local_j < sizeB)
    {
        // Compute global indices
        int global_i = tileA_start + local_i;
        int global_j = tileB_start + local_j;

        // Access the triangles
        TriangleGPU A = trianglesA[local_i];
        TriangleGPU B = trianglesB[local_j];

        // Perform intersection check
        bool intersects = triangleIntersectsGPU(A.A, A.B, A.C, B.A, B.B, B.C);
        if (intersects) {
            // Atomically mark both triangles as intersecting
            atomicExch(&result[global_i], 1);
            atomicExch(&result[global_j], 1);
        }
    }

}


__global__ void checkSingleTriangleKernel(const TriangleGPU T,
                                          const TriangleGPU* trianglesB,
                                          int sizeB,
                                          int T_index,
                                          int* result)
{
    // Calculate the global thread index
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j == T_index)
    {
        // Skip self-intersection
        return;
    }
    
    if (j < sizeB)
    {
        // Access the j-th triangle in B
        TriangleGPU B = trianglesB[j];

        // Perform intersection check between T and B[j]
        bool intersects = triangleIntersectsGPU(T.A, T.B, T.C, B.A, B.B, B.C);

        if (intersects)
        {
            // Atomically mark B[j] as intersecting
            atomicExch(&result[j], 1);
        }
    }

}


__global__ void checkGridBasedKernel(const TriangleGPU* triangles, int numTriangles,
                                     const int* gridCellOffsets, const int* gridCellTriangles,
                                     int3 gridDimensions, double3 minBounds,
                                     double3 cellSize, int* result) {
    int t_id = blockIdx.x * blockDim.x + threadIdx.x;

    if (t_id >= numTriangles) {
        return; // Out-of-bounds threads do nothing
    }

    TriangleGPU T = triangles[t_id];

    double3 tri_min = make_double3(
        fmin(fmin(T.A.x, T.B.x), T.C.x),
        fmin(fmin(T.A.y, T.B.y), T.C.y),
        fmin(fmin(T.A.z, T.B.z), T.C.z)
    );
    double3 tri_max = make_double3(
        fmax(fmax(T.A.x, T.B.x), T.C.x),
        fmax(fmax(T.A.y, T.B.y), T.C.y),
        fmax(fmax(T.A.z, T.B.z), T.C.z)
    );

    int min_cell[3] = {
        max(0, min(gridDimensions.x - 1, int((tri_min.x - minBounds.x) / cellSize.x))),
        max(0, min(gridDimensions.y - 1, int((tri_min.y - minBounds.y) / cellSize.y))),
        max(0, min(gridDimensions.z - 1, int((tri_min.z - minBounds.z) / cellSize.z)))
    };
    int max_cell[3] = {
        max(0, min(gridDimensions.x - 1, int((tri_max.x - minBounds.x) / cellSize.x))),
        max(0, min(gridDimensions.y - 1, int((tri_max.y - minBounds.y) / cellSize.y))),
        max(0, min(gridDimensions.z - 1, int((tri_max.z - minBounds.z) / cellSize.z)))
    };

    // Iterate over grid cells within the bounding box
    for (int x = min_cell[0]; x <= max_cell[0]; x++) {
        for (int y = min_cell[1]; y <= max_cell[1]; y++) {
            for (int z = min_cell[2]; z <= max_cell[2]; z++) {
                // Compute the 1D index for the grid cell
                int cell_index = x + gridDimensions.x * (y + gridDimensions.y * z);

                // Access triangles in this cell using gridCellOffsets
                int start = gridCellOffsets[cell_index];
                int end = gridCellOffsets[cell_index + 1];

                for (int i = start; i < end; i++) {
                    int other_id = gridCellTriangles[i];
                    if (t_id != other_id) {
                        // Access the other triangle
                        TriangleGPU other = triangles[other_id];

                        // Perform triangle intersection check
                        if (triangleIntersectsGPU(T.A, T.B, T.C, other.A, other.B, other.C)) {
                            atomicExch(&result[t_id], 1);
                            atomicExch(&result[other_id], 1);
                        }
                    }
                }
            }
        }
    }
}


