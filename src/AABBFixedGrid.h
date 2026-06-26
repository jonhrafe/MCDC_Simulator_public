#ifndef AABB_FIXED_GRID_H
#define AABB_FIXED_GRID_H


#include <vector>
#include <array>
#include "Eigen/Core"
#include <limits>
#include <cmath>
#include "simerrno.h"
#include "constants.h"   // uint / ulong
#include <iostream>
struct AABB {
    double min_b[3]; // Minimum bounds
    double max_b[3]; // Maximum bounds

    AABB(const Eigen::Vector3d& min_bounds, const Eigen::Vector3d& max_bounds) {
        for (int i = 0; i < 3; ++i) {
            min_b[i] = min_bounds[i];
            max_b[i] = max_bounds[i];
        }
    }

    AABB() {
        for (int i = 0; i < 3; ++i) {
            min_b[i] = std::numeric_limits<double>::max();
            max_b[i] = std::numeric_limits<double>::lowest();
        }
    }
};

class AABBFixedGrid {
public:

    AABBFixedGrid();
    void InitializeGrid(const std::vector<AABB>& aabbs, double cell_size);

    // Method to compute global bounding box for all AABBs
    void computeBoundingBox(const std::vector<AABB>& aabbs);

    // Map a point to a grid cell index
    std::array<int, 3> getCellIndex(const Eigen::Vector3d& point) const;

    // Map an AABB to grid cells
    std::vector<std::array<int, 3>> aabbToGridCells(const AABB& aabb) const;

    // Retrieve AABBs in a specific grid cell
    std::vector<int> getAABBsInCell(const Eigen::Vector3d& point) const;

    std::vector<uint> getAABBsInCells(const AABB& query_aabb) const;

    // Allocation-free variant: fills a caller-owned, reused buffer instead of
    // building an unordered_set + returning a fresh vector on every call (perf B1).
    // No de-duplication: an obstacle spanning several queried cells may appear
    // more than once, which is harmless because the collision handler keeps the
    // nearest hit (idempotent to repeats). 'out' is cleared on entry.
    void getAABBsInCells(const AABB& query_aabb, std::vector<uint>& out) const;

    double computeOptimalCellSize(std::vector<AABB> &aabbs, double memory_limit_mb, double min_cell_size_um) const;

    // Accessors
    std::array<int, 3> getDimensions() const;
    Eigen::Vector3d getMinBounds() const;;
    Eigen::Vector3d getCellSize() const;

    double cell_size;
    std::array<int, 3> grid_dims;
    Eigen::Vector3d min_bounds;
    Eigen::Vector3d max_bounds;

    std::vector<std::vector<int>> grid; // Grid storing AABB indices
    // NOTE: the AABB list is only needed transiently to build the grid; it is NOT
    // stored here (it would be a full second copy -- ~9.3 GB for a 193M-triangle
    // mesh). Queries use 'grid' only.
};

#endif // AABB_FIXED_GRID_H
