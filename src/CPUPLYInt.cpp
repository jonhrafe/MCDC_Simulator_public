#include <iostream>
#include <thread>
#include <fstream>
#include <chrono>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "plyobstacle.h"
#include "FixedGrid.h"

using namespace std;
using namespace Eigen;

void printUsage() {
    cout << "Usage: intersect_detector <input_ply_file> <output_file>" << endl;
    cout << "    <input_ply_file>: Path to the input PLY file containing the mesh." << endl;
    cout << "    <output_file>: Path to the output file where results will be saved." << endl;
}

int main(int argc, char* argv[]) {
    // Check for proper argument count
    if (argc != 3) {
        printUsage();
        return 1;
    }

    // Parse command-line arguments
    string ply_file = argv[1];
    string output_file = argv[2];

    // Load PLY file
    cout << "Loading PLY file: " << ply_file << endl;
    PLYObstacle ply(ply_file);
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
    double min_cell_size = 1.0; // Default minimum cell size
    FixedGrid grid(triangles, min_cell_size);

    cout << "Grid initialized" << endl;
    cout << "Grid size: " << grid.grid_dims[0] << " x " << grid.grid_dims[1] << " x " << grid.grid_dims[2] << endl;

    // Prepare intersection mask
    vector<int> intersected_mask(ply.face_number, 0);

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
            cout << "Nearby triangles: " << nearby_indices.size() << endl;
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
    auto end_time = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(end_time - start_time).count();
    double time_per_triangle = total_time / ply.face_number;

    // Write results to the output file
    ofstream out(output_file);
    if (!out.is_open()) {
        cerr << "Error: Failed to open output file: " << output_file << endl;
        return 1;
    }
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
