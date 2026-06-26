#include "FixedGrid.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>


// Compute the bounding box of all triangles
void FixedGrid::computeBoundingBox() {
    min_bounds = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
    max_bounds = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());

    for (const auto& triangle : triangles) {
        for (const auto& vertex : triangle) {
            min_bounds = min_bounds.cwiseMin(vertex);
            max_bounds = max_bounds.cwiseMax(vertex);
        }
    }
}

// Map a point to a grid cell index
std::array<int, 3> FixedGrid::getCellIndex(const Eigen::Vector3d& point) const {
    std::array<int, 3> cell_index;
    for (int i = 0; i < 3; ++i) {
        cell_index[i] = static_cast<int>((point[i] - min_bounds[i]) / cell_size);
        cell_index[i] = std::clamp(cell_index[i], 0, grid_dims[i] - 1); // Ensure within bounds
    }
    return cell_index;
}

// Find grid cells a triangle overlaps
std::vector<std::array<int, 3>> FixedGrid::triangleToGridCells(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C) const {
    Eigen::Vector3d tri_min = A.cwiseMin(B).cwiseMin(C);
    Eigen::Vector3d tri_max = A.cwiseMax(B).cwiseMax(C);

    // Map bounding box to grid cells
    std::array<int, 3> min_cell = getCellIndex(tri_min);
    std::array<int, 3> max_cell = getCellIndex(tri_max);

    std::vector<std::array<int, 3>> cells;
    for (int x = min_cell[0]; x <= max_cell[0]; ++x) {
        for (int y = min_cell[1]; y <= max_cell[1]; ++y) {
            for (int z = min_cell[2]; z <= max_cell[2]; ++z) {
                cells.push_back({x, y, z});
            }
        }
    }
    return cells;
}

// Constructor
FixedGrid::FixedGrid(const std::vector<std::array<Eigen::Vector3d, 3>>& input_triangles, double min_cell_size)
    : triangles(input_triangles), cell_size(min_cell_size) {
    computeBoundingBox();

    // Compute grid dimensions
    for (int i = 0; i < 3; ++i) {
        grid_dims[i] = static_cast<int>(std::ceil((max_bounds[i] - min_bounds[i]) / cell_size));
    }

    // Initialize grid
    grid.resize(grid_dims[0] * grid_dims[1] * grid_dims[2]);

    // Populate grid
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto& tri = triangles[i];
        auto cells = triangleToGridCells(tri[0], tri[1], tri[2]);
        for (const auto& cell : cells) {
            int index = cell[0] + grid_dims[0] * (cell[1] + grid_dims[1] * cell[2]);
            grid[index].push_back(i);
        }
    }
}

// Get triangles close to a point
std::vector<int> FixedGrid::getTrianglesInCell(const Eigen::Vector3d& point) const {
    auto cell = getCellIndex(point);
    int index = cell[0] + grid_dims[0] * (cell[1] + grid_dims[1] * cell[2]);
    return grid[index];
}


std::vector<int> FixedGrid::getTrianglesForTriangle(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C) const {
    // Compute bounding box of the triangle
    Eigen::Vector3d tri_min = A.cwiseMin(B).cwiseMin(C);
    Eigen::Vector3d tri_max = A.cwiseMax(B).cwiseMax(C);

    // Map bounding box to grid cells
    std::array<int, 3> min_cell = getCellIndex(tri_min);
    std::array<int, 3> max_cell = getCellIndex(tri_max);

    // Collect triangle indices from all relevant cells
    std::unordered_set<int> unique_indices; // Avoid duplicates
    for (int x = min_cell[0]; x <= max_cell[0]; ++x) {
        for (int y = min_cell[1]; y <= max_cell[1]; ++y) {
            for (int z = min_cell[2]; z <= max_cell[2]; ++z) {
                // Compute the 1D index for the grid cell
                int index = x + grid_dims[0] * (y + grid_dims[1] * z);

                // Check if the cell is valid and collect its triangle indices
                if (index >= 0 && index < static_cast<int>(grid.size())) {
                    unique_indices.insert(grid[index].begin(), grid[index].end());
                }
            }
        }
    }

    // Convert the unordered_set to a vector for the return type
    return std::vector<int>(unique_indices.begin(), unique_indices.end());
}


// Convert the grid to CSR (compressed sparse row) format
void FixedGrid::convertGridToCSR() {
    gridOffsets.resize(grid.size() + 1); // One extra for the end offset
    gridTriangles.clear();

    int current_offset = 0;
    for (size_t i = 0; i < grid.size(); ++i) {
        gridOffsets[i] = current_offset;
        gridTriangles.insert(gridTriangles.end(), grid[i].begin(), grid[i].end());
        current_offset += grid[i].size();
    }
    gridOffsets[grid.size()] = current_offset;
}

// CSR (compressed sparse row) grid accessors
size_t FixedGrid::getCellOffsetsSize() const {
    return gridOffsets.size();
}

size_t FixedGrid::getCellTrianglesSize() const {
    return gridTriangles.size();
}

const std::vector<int>& FixedGrid::getCellOffsets() const {
    return gridOffsets;
}

const std::vector<int>& FixedGrid::getCellTriangles() const {
    return gridTriangles;
}

std::array<int, 3> FixedGrid::getDimensions() const {
    return grid_dims;
}

Eigen::Vector3d FixedGrid::getMinBounds() const {
    return min_bounds;
}

Eigen::Vector3d FixedGrid::getCellSize() const {
    return Eigen::Vector3d(cell_size, cell_size, cell_size);
}

