#include <iostream>
#include <thread>
#include <fstream>
#include <chrono>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "dynamicsSimulation.h"
#include "pgsesequence.h"
#include "parallelmcsimulation.h"
#include "voxel.h"
#include "cylinder.h"
#include "simerrno.h"
#include "benchmark.h"
#include "plyobstacle.h"
#include "FixedGrid.h"

typedef unsigned int uint;

using namespace std;
using namespace Eigen;


int main(int argn, char* argv[])
{
    // Load PLY file
    string ply_ = "/home/jonathan/Downloads/output_mesh.ply";
    PLYObstacle ply(ply_);
    cout << "Number of triangles: " << ply.face_number << endl;

    // Create a list of triangles in Eigen format
    vector<array<Vector3d, 3>> triangles;
    for (auto i = 0; i < ply.face_number; i++) {
        Vector3d A, B, C;
        ply.faces[i].getVertex(0, A);
        ply.faces[i].getVertex(1, B);
        ply.faces[i].getVertex(2, C);
        triangles.push_back({A, B, C});
    }

    // Initialize FixedGrid
    double min_cell_size = 1; // Adjust as needed
    FixedGrid grid(triangles, min_cell_size);

    std::cout << "Grid initialized" << std::endl;
    std::cout << "Grid size: " << grid.grid_dims[0] << " x " << grid.grid_dims[1] << " x " << grid.grid_dims[2] << std::endl;

    // Prepare intersection mask
    vector<int> intersected_mask(ply.face_number, 0);

    // Generate output file name with PLY model name appended
    size_t last_slash = ply_.find_last_of('/');
    string ply_name = (last_slash != string::npos) ? ply_.substr(last_slash + 1) : ply_;
    size_t last_dot = ply_name.find_last_of('.');
    ply_name = (last_dot != string::npos) ? ply_name.substr(0, last_dot) : ply_name;

    string output_file = "intersected_faces_CPU_" + ply_name + ".txt";
    ofstream out(output_file);

    // Progress update frequency
    int progress_update = (ply.face_number > 10000) ? 1000 : 10;
    progress_update = (ply.face_number > 100000) ? 10000 : progress_update;

    // Start timing
    auto start_time = chrono::high_resolution_clock::now();

    // Triangle intersection detection using FixedGrid
    for (auto i = 0; i < ply.face_number; i++) {
        // Get triangles near the current triangle
        const auto& tri = triangles[i];
        auto nearby_indices = grid.getTrianglesForTriangle(tri[0], tri[1], tri[2]);

        if ((i % progress_update) == 0) {
            cout << "Progress: " << i << " of " << ply.face_number
                 << " : " << float(i) / float(ply.face_number - 1) * 100 << "%" << endl;

            //sort nearby_indices
            sort(nearby_indices.begin(), nearby_indices.end());

            // for (auto &i : nearby_indices){
            //     cout << i << endl;
            // }
            cout << "Nearby triangles: " << nearby_indices.size() << endl;
            //return 0;
        }

        // Check intersections only with nearby triangles
        for (const auto& j : nearby_indices) {
            if (i != j && ply.faces[i].triangleIntersects(ply.faces[j])) {
                intersected_mask[i] = 1;
                intersected_mask[j] = 1;
            }
        }
    }

    // End timing
    auto end_time     = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(end_time - start_time).count();
    double time_per_triangle = total_time / ply.face_number;

    // Write results to the output file
    for (auto i = 0; i < ply.face_number; i++) {
        out << intersected_mask[i] << endl;
    }
    out.close();

    // Print timing results
    cout << "Results written to " << output_file << endl;
    cout << "Total processing time: " << total_time << " seconds" << endl;
    cout << "Average time per triangle: " << time_per_triangle << " seconds" << endl;

    return 0;
}
