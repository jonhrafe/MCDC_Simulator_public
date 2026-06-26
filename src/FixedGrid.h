#ifndef FIXEDGRID_H
#define FIXEDGRID_H

#include <vector>
#include <array>
#include <Eigen/Dense> // For Eigen::Vector3d

class FixedGrid {
public:
    std::vector<std::array<Eigen::Vector3d, 3>> triangles; // List of triangles (A, B, C)
    std::vector<std::vector<int>> grid;                   // Grid storing triangle indices
    std::vector<int> gridOffsets;                         // Start index of triangles for each cell (CSR format)
    std::vector<int> gridTriangles;                       // Flat array of triangle indices (CSR format)

    double cell_size;                                     // Minimum cell size

    // Compute the bounding box of all triangles
    void computeBoundingBox();

    // Map a point to a grid cell index
    std::array<int, 3> getCellIndex(const Eigen::Vector3d& point) const;

    // Find grid cells a triangle overlaps
    std::vector<std::array<int, 3>> triangleToGridCells(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C) const;

    std::array<int, 3> grid_dims;                         // Number of cells in each dimension
    Eigen::Vector3d min_bounds;                           // Bounding box minimum (x, y, z)  
    Eigen::Vector3d max_bounds;                           // Bounding box maximum (x, y, z)

    // Constructor: initializes the grid with triangles and minimum cell size
    FixedGrid(const std::vector<std::array<Eigen::Vector3d, 3>>& input_triangles, double min_cell_size);

    // Get triangles close to a given point
    std::vector<int> getTrianglesInCell(const Eigen::Vector3d& point) const;

    std::vector<int> getTrianglesForTriangle(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C) const;

    // Flat CSR (compressed sparse row) grid accessors
    size_t getCellOffsetsSize() const;                   // Returns size of gridOffsets
    size_t getCellTrianglesSize() const;                 // Returns size of gridTriangles
    const std::vector<int>& getCellOffsets() const;      // Returns the gridOffsets array
    const std::vector<int>& getCellTriangles() const;    // Returns the gridTriangles array
    std::array<int, 3> getDimensions() const;            // Returns grid dimensions
    Eigen::Vector3d getMinBounds() const;                // Returns grid minimum bounds
    Eigen::Vector3d getCellSize() const;                 // Returns cell size

    void convertGridToCSR();
};

#endif // FIXEDGRID_H
