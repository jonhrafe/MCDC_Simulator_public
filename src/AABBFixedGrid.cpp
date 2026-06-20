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
    // Do NOT copy 'aabbs' into a member: it is only needed to build the grid here,
    // and a second resident copy is ~9.3 GB for a 193M-triangle mesh.
    this->cell_size = cell_size;

    if (cell_size <= 0) {
        throw std::invalid_argument("Cell size must be positive.");
    }

    computeBoundingBox(aabbs);

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

void AABBFixedGrid::computeBoundingBox(const std::vector<AABB>& aabbs) {
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

        // Check if the cell is completely outside the grid
    for (int i = 0; i < 3; ++i) {
        if (cell[i] < 0 || cell[i] >= grid_dims[i]) {
            return {}; // Return empty if point is out of bounds
        }
    }
    return grid[index];
}

void AABBFixedGrid::getAABBsInCells(const AABB& query_aabb, std::vector<uint>& out) const {
    out.clear();
    std::array<int, 3> min_cell = getCellIndex(Eigen::Vector3d(query_aabb.min_b[0], query_aabb.min_b[1], query_aabb.min_b[2]));
    std::array<int, 3> max_cell = getCellIndex(Eigen::Vector3d(query_aabb.max_b[0], query_aabb.max_b[1], query_aabb.max_b[2]));

    for (int x = min_cell[0]; x <= max_cell[0]; ++x) {
        for (int y = min_cell[1]; y <= max_cell[1]; ++y) {
            for (int z = min_cell[2]; z <= max_cell[2]; ++z) {
                int index = x + grid_dims[0] * (y + grid_dims[1] * z);
                if (index >= 0 && index < static_cast<int>(grid.size())) {
                    const std::vector<int>& cell = grid[index];
                    out.insert(out.end(), cell.begin(), cell.end());
                }
            }
        }
    }
}

std::vector<uint> AABBFixedGrid::getAABBsInCells(const AABB& query_aabb) const {
    // Back-compat wrapper (de-duplicated). The per-step hot path uses the buffer
    // overload above; this is kept for any non-hot caller.
    std::vector<uint> tmp;
    getAABBsInCells(query_aabb, tmp);
    std::unordered_set<uint> unique_indices(tmp.begin(), tmp.end());
    return std::vector<uint>(unique_indices.begin(), unique_indices.end());
}

double AABBFixedGrid::computeOptimalCellSize(std::vector<AABB> &aabbs, double memory_limit_mb, double min_cell_size_um) const
{
    
    if(aabbs.empty()){
        // return a big number so the cells are not created
        return 1000000;
    }   

    // Convert memory limit to bytes
    const double memory_limit_bytes = memory_limit_mb * 1024 * 1024;

    // Compute global bounding box volume
    Eigen::Vector3d global_min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
    Eigen::Vector3d global_max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());

    for (const auto& aabb : aabbs) {
        for (int i = 0; i < 3; ++i) {
            global_min[i] = std::min(global_min[i], aabb.min_b[i]);
            global_max[i] = std::max(global_max[i], aabb.max_b[i]);
        }
    }
    Eigen::Vector3d global_size = global_max - global_min;
    double global_volume = global_size.prod();

    // Compute average AABB volume
    double total_aabb_volume = 0.0;
    for (const auto& aabb : aabbs) {
        double aabb_volume = 1.0;
        for (int i = 0; i < 3; ++i) {
            aabb_volume *= (aabb.max_b[i] - aabb.min_b[i]);
        }
        total_aabb_volume += aabb_volume;
    }
    double avg_aabb_volume = total_aabb_volume / aabbs.size();

    // Estimate grid dimensions based on global volume and target memory limit
    double estimated_cell_size = std::cbrt(avg_aabb_volume)/2; // 1 AABBs per cell heuristic

    // Adjust cell size to fit within memory constraints
    double grid_cells = (global_volume / std::pow(estimated_cell_size, 3));
    double estimated_memory_usage = grid_cells * sizeof(std::vector<int>); // Approx memory per cell

    while (estimated_memory_usage > memory_limit_bytes && estimated_cell_size < global_size.minCoeff()) {
        estimated_cell_size *= 1.05; // Gradually increase cell size to reduce memory usage
        grid_cells = (global_volume / std::pow(estimated_cell_size, 3));
        estimated_memory_usage = grid_cells * sizeof(std::vector<int>);

        if(estimated_memory_usage > memory_limit_bytes){
            std::string message = "AABB collision max memory reached. (" + std::to_string(memory_limit_mb) + " MB) increase mem limit in constants.h \n";
            SimErrno::info(message,std::cout);
        }
    }

    // Enforce minimum cell size
    double min_cell_size_mm = min_cell_size_um * 1e-3; // Convert microns to mm
    return std::max(estimated_cell_size, min_cell_size_mm);
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
