#include "AABBFixedGrid.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

AABBFixedGrid::AABBFixedGrid() {
    cell_size = 0;
    grid_dims.fill(0);
    min_bounds = Eigen::Vector3d::Zero();
    max_bounds = Eigen::Vector3d::Zero();
}


void AABBFixedGrid::InitializeGrid(const std::vector<AABB>& aabbs, double cell_size){
    this->aabbs = aabbs;
    this->cell_size = cell_size;

    if (cell_size <= 0) {
        throw std::invalid_argument("Cell size must be positive.");
    }

    computeBoundingBox();

    // Compute grid dimensions
    for (int i = 0; i < 3; ++i) {
        grid_dims[i] = static_cast<int>(std::ceil((max_bounds[i] - min_bounds[i]) / cell_size));
    }

    // Initialize grid
    grid.resize(grid_dims[0] * grid_dims[1] * grid_dims[2]);

    // Populate grid with AABBs
    for (size_t i = 0; i < aabbs.size(); ++i) {
        auto cells = aabbToGridCells(aabbs[i]);
        for (const auto& cell : cells) {
            int index = cell[0] + grid_dims[0] * (cell[1] + grid_dims[1] * cell[2]);
            grid[index].push_back(i);
        }
    }
}

void AABBFixedGrid::computeBoundingBox() {
    min_bounds = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
    max_bounds = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());

    for (const auto& aabb : aabbs) {
        for (int i = 0; i < 3; ++i) {
            min_bounds[i] = std::min(min_bounds[i], aabb.min_b[i]);
            max_bounds[i] = std::max(max_bounds[i], aabb.max_b[i]);
        }
    }
}

std::array<int, 3> AABBFixedGrid::getCellIndex(const Eigen::Vector3d& point) const {
    std::array<int, 3> cell_index;
    for (int i = 0; i < 3; ++i) {
        cell_index[i] = static_cast<int>((point[i] - min_bounds[i]) / cell_size);
        cell_index[i] = std::clamp(cell_index[i], 0, grid_dims[i] - 1); // Ensure within bounds
    }
    return cell_index;
}

std::vector<std::array<int, 3>> AABBFixedGrid::aabbToGridCells(const AABB& aabb) const {
    std::array<int, 3> min_cell = getCellIndex(Eigen::Vector3d(aabb.min_b[0], aabb.min_b[1], aabb.min_b[2]));
    std::array<int, 3> max_cell = getCellIndex(Eigen::Vector3d(aabb.max_b[0], aabb.max_b[1], aabb.max_b[2]));

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

std::vector<int> AABBFixedGrid::getAABBsInCell(const Eigen::Vector3d& point) const {
    auto cell = getCellIndex(point);
    int index = cell[0] + grid_dims[0] * (cell[1] + grid_dims[1] * cell[2]);
    return grid[index];
}

std::vector<uint> AABBFixedGrid::getAABBsInCells(const AABB& query_aabb) const {
    std::array<int, 3> min_cell = getCellIndex(Eigen::Vector3d(query_aabb.min_b[0], query_aabb.min_b[1], query_aabb.min_b[2]));
    std::array<int, 3> max_cell = getCellIndex(Eigen::Vector3d(query_aabb.max_b[0], query_aabb.max_b[1], query_aabb.max_b[2]));

    std::unordered_set<uint> unique_indices;

    for (int x = min_cell[0]; x <= max_cell[0]; ++x) {
        for (int y = min_cell[1]; y <= max_cell[1]; ++y) {
            for (int z = min_cell[2]; z <= max_cell[2]; ++z) {
                int index = x + grid_dims[0] * (y + grid_dims[1] * z);
                if (index >= 0 && index < static_cast<int>(grid.size())) {
                    unique_indices.insert(grid[index].begin(), grid[index].end());
                }
            }
        }
    }

    return std::vector<uint>(unique_indices.begin(), unique_indices.end());
}

std::array<int, 3> AABBFixedGrid::getDimensions() const {
    return grid_dims;
}

Eigen::Vector3d AABBFixedGrid::getMinBounds() const {
    return min_bounds;
}

Eigen::Vector3d AABBFixedGrid::getCellSize() const {
    return Eigen::Vector3d(cell_size, cell_size, cell_size);
}
