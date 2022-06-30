#pragma once

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

using namespace std;

// Specifies the key of containers (or bins in the grid)
typedef size_t ContainerKey;

// Specifies the identification number of cells (must be sequential from 0 to ncell - 1)
typedef size_t CellId;

// Defines the container type
typedef set<CellId> Container_t;

// Defines the containers type
typedef map<ContainerKey, Container_t> Containers_t;

// Defines the bounding box of a cell
const size_t N_MIN_MAX = 2;  // 2 means {min,max}
const size_t I_MIN = 0;
const size_t I_MAX = 1;
typedef vector<vector<double>> BboxMinMax;  // [ndim][N_MIN_MAX]

// Defines a tool to search the cell where a point is located within a mesh
//
// # Reference
//
// * Durand, Farias, and Pedroso (2015) Computing intersections between
//   non-compatible curves and finite elements, Computational Mechanics;
//   DOI=10.1007/s00466-015-1181-y
struct GridSearch {
    size_t ndim;                        // space dimension
    vector<size_t> ndiv;                // (ndim) number of divisions along each direction
    vector<double> xmin;                // (ndim) min values
    vector<double> xmax;                // (ndim) max values
    double side_length;                 // (ndim) side lengths of each container
    vector<BboxMinMax> bounding_boxes;  // (ncell) bounding boxes
    Containers_t containers;            // holds all items

    // Calculates the key of the container where the point should fall in
    static inline ContainerKey calc_container_key(
        size_t ndim,
        double side_length,
        vector<size_t> const &ndiv,
        vector<double> const &xmin,
        vector<double> const &x) {
        size_t ix = static_cast<size_t>((x[0] - xmin[0]) / side_length);  // (Eq. 8)
        size_t iy = static_cast<size_t>((x[1] - xmin[1]) / side_length);
        if (ix == ndiv[0]) {
            ix -= 1;  // point is on max edge => move to inner container
        }
        if (iy == ndiv[1]) {
            iy -= 1;  // point is on max edge => move to inner container
        }
        if (ndim == 2) {
            return ix + iy * ndiv[0];
        }
        size_t iz = static_cast<size_t>((x[2] - xmin[2]) / side_length);
        if (iz == ndiv[2]) {
            iz -= 1;  // point is on max edge => move to inner container
        }
        return ix + iy * ndiv[0] + iz * ndiv[0] * ndiv[1];
    }

    // Allocates new instance
    ///
    /// # Input
    ///
    /// * `ndim` -- is the space dimension (2 or 3)
    /// * `ncell` -- is the number of cells (e.g., triangle/tetrahedron) in the mesh.
    ///     - All cells are numbered from `0` to `ncell - 1`
    ///     - The index of the cell in a mesh is also called CellId (`cell_id`)
    /// * `get_nnode` -- is a function of the `cell_id` that returns the number of nodes `nnode` of the cell
    /// * `get_x` -- is a function of the `cell_id` and the local index of the node/point `m`.
    ///    This function returns the coordinates `x` of the point.
    /// * `tolerance` -- is a tolerance to expand the bounding box of cells and compare points; e.g. 1e-4
    ///     - If None, [GS_DEFAULT_TOLERANCE] is used
    /// * `border_tol` -- is a tolerance used to expand the border a little bit and then
    ///   accommodate eventual imprecision near the borders; e.g. 1e-2
    ///     - If None, [GS_DEFAULT_BORDER_TOL] is used
    static std::unique_ptr<GridSearch> make_new(
        size_t ndim,
        size_t ncell,
        function<size_t(size_t)> get_nnode,
        function<vector<double> const *(size_t, size_t)> get_x);
    // size_t (*get_nnode)(size_t),
    // vector<double> const &(*get_x)(size_t, size_t));

    // Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
    //
    // * Returns the CellId or -1 if no cell contains the point
    // * `is_in_cell` is a function of cell_id and x that tells whether the point si in the cell or not
    int find_cell(vector<double> &x, function<bool(size_t, vector<double> const *)> is_in_cell);

    // Print details about the grid (e.g., for debugging)
    void print_details();
};