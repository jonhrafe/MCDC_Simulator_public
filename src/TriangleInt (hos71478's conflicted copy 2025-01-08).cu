// TriangleInt.cu

#include "CUDAInt.cuh"
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include "triangle.h"
#include "plyobstacle.h"
#include <Eigen/Dense>  // Ensure Eigen is installed and included properly
#include <chrono>
#include "FixedGrid.h"

// Host-side CUDA error checking macros
#define CUDA_ERROR_CHECK

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }

    // More careful checking. However, this will affect performance.
    // Comment away if needed.
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
    return;
}

// Utility to cast Eigen::Vector3d to double3 for CUDA
double3 toDouble3(const Eigen::Vector3d& v) {
    return {v.x(), v.y(), v.z()};
}


using namespace std;
using namespace Eigen;

int main(int argn, char* argv[])
{
    // Load PLY file
    string ply_ = "/home/jonathan/Downloads/output_mesh.ply"; //hp_outer.ply"; // output_mesh.ply
    PLYObstacle ply(ply_);
    cout << "Number of triangles: " << ply.face_number << endl;

    // Create a list of triangles in Eigen format
    vector<array<Vector3d, 3>> triangles;
    vector<TriangleGPU> h_triangles; // For CUDA
    for (auto i = 0; i < ply.face_number; i++) {
        Vector3d A, B, C;
        ply.faces[i].getVertex(0, A);
        ply.faces[i].getVertex(1, B);
        ply.faces[i].getVertex(2, C);
        triangles.push_back({A, B, C});

        // Add to CUDA-friendly structure
        TriangleGPU t;
        t.A = {A.x(), A.y(), A.z()};
        t.B = {B.x(), B.y(), B.z()};
        t.C = {C.x(), C.y(), C.z()};
        h_triangles.push_back(t);
    }

    // Initialize FixedGrid
    double min_cell_size = 1.0; // Adjust as needed
    FixedGrid grid(triangles, min_cell_size);

    // Prepare grid for CUDA
    grid.convertGridToCSR();

    // Output file setup
    size_t last_slash = ply_.find_last_of('/');
    string ply_name = (last_slash != string::npos) ? ply_.substr(last_slash + 1) : ply_;
    size_t last_dot = ply_name.find_last_of('.');
    ply_name = (last_dot != string::npos) ? ply_name.substr(0, last_dot) : ply_name;
    string output_file = "intersected_faces_GPU_" + ply_name + ".txt";
    ofstream out(output_file);

    // Prepare results array
    vector<int> h_result(ply.face_number, 0); // Host result
    int* d_result;
    CudaSafeCall(cudaMalloc(&d_result, h_result.size() * sizeof(int)));
    CudaSafeCall(cudaMemset(d_result, 0, h_result.size() * sizeof(int)));

    // Allocate and copy grid data to device
    int* d_gridCellOffsets;
    int* d_gridCellTriangles;
    CudaSafeCall(cudaMalloc(&d_gridCellOffsets, grid.getCellOffsetsSize() * sizeof(int)));
    CudaSafeCall(cudaMalloc(&d_gridCellTriangles, grid.getCellTrianglesSize() * sizeof(int)));
    CudaSafeCall(cudaMemcpy(d_gridCellOffsets, grid.getCellOffsets().data(), grid.getCellOffsetsSize() * sizeof(int), cudaMemcpyHostToDevice));
    CudaSafeCall(cudaMemcpy(d_gridCellTriangles, grid.getCellTriangles().data(), grid.getCellTrianglesSize() * sizeof(int), cudaMemcpyHostToDevice));


    // Allocate and copy triangle data to device
    TriangleGPU* d_triangles;
    cudaMalloc(&d_triangles, h_triangles.size() * sizeof(TriangleGPU));
    cudaMemcpy(d_triangles, h_triangles.data(), h_triangles.size() * sizeof(TriangleGPU), cudaMemcpyHostToDevice);

    // Prepare grid dimensions and bounds

    double3 minBounds   = make_double3(grid.getMinBounds().x(), grid.getMinBounds().y(), grid.getMinBounds().z());
    double3 cellSize    = make_double3(grid.getCellSize().x(), grid.getCellSize().y(), grid.getCellSize().z());
    int3 gridDimensions = make_int3(grid.getDimensions()[0], grid.getDimensions()[1], grid.getDimensions()[2]);

    // Start timing
    auto start_time = chrono::high_resolution_clock::now();

    // Configure kernel launch parameters
    int threadsPerBlock = 256;
    int numBlocks = (ply.face_number + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the kernel
    checkGridBasedKernel<<<numBlocks, threadsPerBlock>>>(
        d_triangles, ply.face_number,
        d_gridCellOffsets, d_gridCellTriangles,
        gridDimensions, minBounds, cellSize,
        d_result
    );

    cudaDeviceSynchronize();
    // Add CUDA-safe calls to check for errors
    CudaCheckError();


    // End timing
    auto end_time     = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(end_time - start_time).count();
    double time_per_triangle = total_time / ply.face_number;

    // Copy results back to host
    cudaMemcpy(h_result.data(), d_result, h_result.size() * sizeof(int), cudaMemcpyDeviceToHost);

    // Write results to the output file
    for (auto i = 0; i < ply.face_number; i++) {
        out << h_result[i] << endl;
    }
    out.close();

    // Print timing results
    cout << "Results written to " << output_file << endl;
    cout << "Total processing time: " << total_time << " seconds" << endl;
    cout << "Average time per triangle: " << time_per_triangle << " seconds" << endl;

    // Free device memory
    cudaFree(d_result);
    cudaFree(d_gridCellOffsets);
    cudaFree(d_gridCellTriangles);
    cudaFree(d_triangles);

    return 0;
}

//     // 3) Prepare result array
//     std::vector<int> h_result(N, 0); // Initialize with zeros

//     // 4) Allocate device memory and copy data
//     TriangleGPU* d_triangles = nullptr;
//     int* d_result = nullptr;

//     // Allocate device memory
//     CudaSafeCall(cudaMalloc(&d_triangles, N * sizeof(TriangleGPU)));
//     CudaSafeCall(cudaMalloc(&d_result, N * sizeof(int)));

//     // Copy triangles to device
//     CudaSafeCall(cudaMemcpy(d_triangles, h_triangles.data(), N * sizeof(TriangleGPU), cudaMemcpyHostToDevice));

//     // Initialize result array on device to zero
//     CudaSafeCall(cudaMemset(d_result, 0, N * sizeof(int)));

//     // 5) Define grid and block dimensions

//     const int TILE_SIZE = 100000; // Adjust based on GPU capabilities
//     int num_tiles = (N + TILE_SIZE - 1) / TILE_SIZE;

//     dim3 block(16, 16);
//     dim3 grid( (TILE_SIZE + block.x - 1)/block.x,
//                (TILE_SIZE + block.y - 1)/block.y );

//     std::cout << "Launching kernel " << std::endl;

//     // <-- Added: Create CUDA events for timing
//     cudaEvent_t start, stop;
//     CudaSafeCall(cudaEventCreate(&start));
//     CudaSafeCall(cudaEventCreate(&stop));

//     // <-- Added: Record the start event
//     CudaSafeCall(cudaEventRecord(start, 0));

//     // 6) Launch the kernel
//     // Iterate over all tile pairs
//     for(int tileA = 0; tileA < num_tiles; ++tileA)
//     {
//         int tileA_start = tileA * TILE_SIZE;
//         int sizeA = min(TILE_SIZE, N - tileA_start);

//         std::cout << "Processing tile " << tileA+1 << " of " << num_tiles << std::endl;

//         for(int tileB = tileA; tileB < num_tiles; ++tileB)
//         {
//             int tileB_start = tileB * TILE_SIZE;
//             int sizeB = min(TILE_SIZE, N - tileB_start);

//             // Define grid dimensions based on current tile sizes
//             dim3 current_grid( (sizeA + block.x - 1) / block.x,
//                             (sizeB + block.y - 1) / block.y );

//             // Pointers to the current tile subsets
//             const TriangleGPU* current_trianglesA = d_triangles + tileA_start;
//             const TriangleGPU* current_trianglesB = d_triangles + tileB_start;

//             // Launch the kernel
//             checkAllPairsKernel<<<current_grid, block>>>(current_trianglesA, sizeA,
//                                                         current_trianglesB, sizeB,
//                                                         tileA_start, tileB_start,
//                                                         d_result);
//             // Check for kernel launch errors
//             CudaCheckError();
//         }
//     }
//     // <-- Added: Record the stop event
//     CudaSafeCall(cudaEventRecord(stop, 0));

//     // <-- Added: Wait for the stop event to complete
//     CudaSafeCall(cudaEventSynchronize(stop));

//     // <-- Added: Calculate elapsed time
//     float milliseconds = 0;
//     CudaSafeCall(cudaEventElapsedTime(&milliseconds, start, stop));

//     // <-- Added: Destroy CUDA events
//     CudaSafeCall(cudaEventDestroy(start));
//     CudaSafeCall(cudaEventDestroy(stop));

//     // Compute triangles per second
//     // Interpretation: Number of triangle pairs processed per second
//     // Total triangle pairs = N * N
//     // Alternatively, if you consider each triangle being processed individually, adjust accordingly
//     double time_seconds = milliseconds / 1000.0;
//     double triangle_pairs = static_cast<double>(N) * static_cast<double>(N);
//     double triangles_per_second = triangle_pairs / time_seconds;

//     std::cout << "Kernel execution time: " << milliseconds << " ms" << std::endl;
//     std::cout << "Triangles processed per second: " << triangles_per_second << std::endl;

//     // 7) Copy results back to host
//     CudaSafeCall(cudaMemcpy(h_result.data(), d_result, N * sizeof(int), cudaMemcpyDeviceToHost));

//     // 8) Write results to a file
//     std::string output_file = "intersected_faces_cuda.txt";
//     std::ofstream out(output_file);
//     if (!out.is_open()) {
//         std::cerr << "Failed to open output file: " << output_file << std::endl;
//         // Free device memory before exiting
//         cudaFree(d_triangles);
//         cudaFree(d_result);
//         return -1;
//     }

//     for (int i = 0; i < N; i++) {
//         out << h_result[i] << std::endl;
//     }

//     out.close();
//     std::cout << "Results written to " << output_file << std::endl;

//     // 9) Free device memory
//     CudaSafeCall(cudaFree(d_triangles));
//     CudaSafeCall(cudaFree(d_result));

//     return 0;
// }
