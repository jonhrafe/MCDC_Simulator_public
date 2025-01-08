#include "CUDAInt.cuh"
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include "triangle.h"
#include "plyobstacle.h"
#include "Eigen/Dense" 
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

void printUsage() {
    std::cout << "Usage: intersect_detector <input_ply_file> <output_file>" << std::endl;
    std::cout << "    <input_ply_file>: Path to the input PLY file containing the mesh." << std::endl;
    std::cout << "    <output_file>: Path to the output file where results will be saved." << std::endl;
}

int main(int argc, char* argv[]) {
    // Check for proper argument count
    if (argc != 3) {
        printUsage();
        return 1;
    }

    // Parse command-line arguments
    std::string ply_file = argv[1];
    std::string output_file = argv[2];

    // Load PLY file
    std::cout << "Loading PLY file: " << ply_file << std::endl;
    PLYObstacle ply(ply_file);
    std::cout << "Number of triangles: " << ply.face_number << std::endl;

    // Create a list of triangles in Eigen format
    std::vector<std::array<Eigen::Vector3d, 3>> triangles;
    std::vector<TriangleGPU> h_triangles; // For CUDA
    for (auto i = 0; i < ply.face_number; i++) {
        Eigen::Vector3d A, B, C;
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
    grid.convertGridToCSR();

    // Prepare results array
    std::vector<int> h_result(ply.face_number, 0); // Host result
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
    CudaSafeCall(cudaMalloc(&d_triangles, h_triangles.size() * sizeof(TriangleGPU)));
    CudaSafeCall(cudaMemcpy(d_triangles, h_triangles.data(), h_triangles.size() * sizeof(TriangleGPU), cudaMemcpyHostToDevice));

    // Prepare grid dimensions and bounds
    double3 minBounds = make_double3(grid.getMinBounds().x(), grid.getMinBounds().y(), grid.getMinBounds().z());
    double3 cellSize = make_double3(grid.getCellSize().x(), grid.getCellSize().y(), grid.getCellSize().z());
    int3 gridDimensions = make_int3(grid.getDimensions()[0], grid.getDimensions()[1], grid.getDimensions()[2]);

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Configure kernel launch parameters
    int threadsPerBlock = 128;
    int numBlocks = (ply.face_number + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the kernel
    checkGridBasedKernel<<<numBlocks, threadsPerBlock>>>(
        d_triangles, ply.face_number,
        d_gridCellOffsets, d_gridCellTriangles,
        gridDimensions, minBounds, cellSize,
        d_result
    );
    CudaCheckError();
    cudaDeviceSynchronize();
    CudaCheckError();

    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(end_time - start_time).count();
    double time_per_triangle = total_time / ply.face_number;

    // Copy results back to host
    CudaSafeCall(cudaMemcpy(h_result.data(), d_result, h_result.size() * sizeof(int), cudaMemcpyDeviceToHost));

    // Write results to the output file
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Failed to open output file: " << output_file << std::endl;
        return 1;
    }
    for (auto i = 0; i < ply.face_number; i++) {
        out << h_result[i] << std::endl;
    }
    out.close();

    // Print timing results
    std::cout << "Results written to " << output_file << std::endl;
    std::cout << "Total processing time: " << total_time << " seconds" << std::endl;
    std::cout << "Average time per triangle: " << time_per_triangle << " seconds" << std::endl;

    // Free device memory
    cudaFree(d_result);
    cudaFree(d_gridCellOffsets);
    cudaFree(d_gridCellTriangles);
    cudaFree(d_triangles);

    return 0;
}
