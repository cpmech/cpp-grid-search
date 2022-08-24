#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "tritet_coords.h"

using namespace std;

// Indicates that the find function found a cell with success
const size_t SUCCESS = 1000;

// Indicates that the find function could not find a cell
const size_t CELL_NOT_FOUND = 2000;

// Indicates that the find function failed because the point is outside the grid
const size_t FAIL_POINT_OUTSIDE = 3000;

// GridSearch tolerance for all directions
const double GS_TOLERANCE = 1e-4;

// GridSearch border tolerance to handle imprecision near the borders
const double GS_BORDER_TOL = 1e-2;

/// Defines a bounding-box ratio coefficient in [0,1] to mark cells as "large" cells compared the others
///
/// Example: If a cell has a max(bounding box length) / largest >= 0.75,
///          Then this cell is put into the large_cells list
///
/// This constant enables the strategy of putting some cells in a selected list of large cells.
/// The other constant [GS_LARGE_CELLS_MAX_COUNT_PCT] also influences the algorithm by enabling or disabling this strategy.
const double GS_LARGE_CELLS_CUTOFF = 0.75;

/// Controls the count percentage in [0,1] of large cells allowed in the selected array of large cells.
///
/// Example: If 20% of cells are too large, put them in a separated list of large cells;
///          otherwise, the mesh is homogeneous enough so put all cells in the containers.
///
/// If the actual ratio is greater than this constant (i.e., many cells are large such as in a homogeneous mesh),
/// the "strategy" of selecting large cells is abandoned. The other constant [GS_LARGE_CELLS_CUTOFF] also influences the algorithm.
const double GS_LARGE_CELLS_MAX_COUNT_PCT = 0.2;

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

// Defines a tool to search the cell (triangle or tetrahedron) where a point is located within a mesh
//
// # Reference
//
// * Durand, Farias, and Pedroso (2015) Computing intersections between
//   non-compatible curves and finite elements, Computational Mechanics;
//   DOI=10.1007/s00466-015-1181-y
struct GridSearch {
    size_t ndim;                           // space dimension
    vector<vector<double>> const &points;  // holds a reference to the list of points
    vector<vector<size_t>> const &cells;   // holds a reference to the list of cells
    vector<size_t> ndiv;                   // (ndim) number of divisions along each direction
    vector<double> xmin;                   // (ndim) min values
    vector<double> xmax;                   // (ndim) max values
    double side_length;                    // side length of a container
    vector<BboxMinMax> bounding_boxes;     // (ncell) bounding boxes
    Containers_t containers;               // holds all items
    Container_t large_cells;               // holds the ids of large cells

    // Calculates the key of the container where the point should fall in
    static inline ContainerKey calc_container_key(size_t ndim,
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
    //
    // # Input
    //
    // * `points` -- is a list of coordinates such as
    //   `[[x0,y0,T0], [x1,y1,T1], [x2,y2,T2], [x3,y3,T3]]` in 2D, or
    //   `[[x0,y0,z0,T0], [x1,y1,z1,T1], [x2,y2,z2,T2], [x3,y3,z3,T3]]` in 3D.
    //    The last value `Ti` is optional.
    // * `cells` -- is a list of connectivity (topology) such as
    //   `[[0,2,1], [2,1,0]]` in 2D, or
    //   `[[0,2,1,3], [2,1,0,3]]` in 3D
    static std::unique_ptr<GridSearch> make_new(size_t ndim,
                                                vector<vector<double>> const &points,
                                                vector<vector<size_t>> const &cells) {
        // check
        if (ndim < 2 || ndim > 3) {
            throw "ndim must be 2 or 3";
        }
        size_t npoint = points.size();
        if (npoint < ndim + 1) {
            throw "number of points must be >= ndim + 1";
        }
        size_t ncol = points[0].size();
        if (ncol != ndim + 1) {
            throw "points.ncol must be = ndim + 1";
        }
        size_t ncell = cells.size();
        if (ncell < 1) {
            throw "number of cells must be >= 1";
        }
        size_t nnode = cells[0].size();
        if (nnode != ndim + 1) {
            throw "number of cell nodes (nnode) must be = ndim + 1";
        }

        // allocate some variables
        auto x = vector<double>(ndim);
        auto xmin = vector<double>(ndim);
        auto xmax = vector<double>(ndim);
        auto x_min_max = BboxMinMax(ndim);
        auto bbox_largest = vector<double>(ndim);
        vector<BboxMinMax> bounding_boxes;
        for (size_t i = 0; i < ndim; i++) {
            xmin[i] = numeric_limits<double>::max();
            xmax[i] = numeric_limits<double>::min();
            x_min_max[i] = vector<double>(N_MIN_MAX);
            bbox_largest[i] = numeric_limits<double>::min();
        }

        // find limits, bounding boxes, and largest bounding box
        for (size_t id = 0; id < ncell; id++) {
            for (size_t m = 0; m < nnode; m++) {
                size_t p = cells[id][m];
                for (size_t i = 0; i < ndim; i++) {
                    // coords
                    x[i] = points[p][i];
                    // limits
                    xmin[i] = min(xmin[i], x[i]);
                    xmax[i] = max(xmax[i], x[i]);
                    // bounding box
                    if (m == 0) {
                        x_min_max[i][I_MIN] = x[i];
                        x_min_max[i][I_MAX] = x[i];
                    } else {
                        x_min_max[i][I_MIN] = min(x_min_max[i][I_MIN], x[i]);
                        x_min_max[i][I_MAX] = max(x_min_max[i][I_MAX], x[i]);
                    }
                }
            }
            // largest bounding box
            for (size_t i = 0; i < ndim; i++) {
                bbox_largest[i] = max(bbox_largest[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            // add to bounding box maps
            bounding_boxes.push_back(x_min_max);
        }

        // compute the largest length of largest bounding box
        double bbox_largest_length = bbox_largest[0];
        for (size_t i = 1; i < ndim; i++) {
            bbox_largest_length = max(bbox_largest_length, bbox_largest[i]);
        }

        // handle large cells
        double bbox_large_length = numeric_limits<double>::min();  // largest length lower than the cutoff
        Container_t large_cells;
        for (size_t id = 0; id < ncell; id++) {
            auto x_min_max = bounding_boxes[id];
            double largest_length = x_min_max[0][I_MAX] - x_min_max[0][I_MIN];
            for (size_t i = 1; i < ndim; i++) {
                largest_length = max(largest_length, x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            if (largest_length >= GS_LARGE_CELLS_CUTOFF * bbox_largest_length) {
                large_cells.insert(id);
            } else {
                for (size_t i = 0; i < ndim; i++) {
                    bbox_large_length = max(bbox_large_length, x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
                }
            }
        }

        // make the side_length equal to the largest bounding box dimension
        double pct = static_cast<double>(large_cells.size()) / static_cast<double>(ncell);
        double side_length;
        if (pct <= GS_LARGE_CELLS_MAX_COUNT_PCT) {
            side_length = bbox_large_length;  // use the largest length lower than the cutoff
        } else {
            large_cells.clear();                // abandon the large cells strategy
            side_length = bbox_largest_length;  // use the largest length among them all
        };

        // expand side_length by two times the tolerance (left and right)
        side_length += 2.0 * GS_TOLERANCE;

        // expand borders
        for (size_t i = 0; i < ndim; i++) {
            xmin[i] -= GS_BORDER_TOL;
            xmax[i] += GS_BORDER_TOL;
        }

        // number of divisions
        vector<size_t> ndiv(ndim);
        for (size_t i = 0; i < ndim; i++) {
            ndiv[i] = ceil((xmax[i] - xmin[i]) / side_length);
        }

        // update xmax after deciding on the side_length and number of divisions
        for (size_t i = 0; i < ndim; i++) {
            xmax[i] = xmin[i] + side_length * ((double)ndiv[i]);
        }

        // insert cells
        Containers_t containers;
        for (size_t id = 0; id < ncell; id++) {
            auto iter = large_cells.find(id);
            if (iter != large_cells.end()) {
                continue;  // skip large cells
            }
            auto x_min_max = bounding_boxes[id];
            for (size_t r = 0; r < N_MIN_MAX; r++) {
                for (size_t s = 0; s < N_MIN_MAX; s++) {
                    for (size_t t = 0; t < ndim - 1; t++) {
                        x[0] = x_min_max[0][r];
                        x[1] = x_min_max[1][s];
                        if (ndim == 3) {
                            x[2] = x_min_max[2][t];
                        }
                        auto key = calc_container_key(ndim, side_length, ndiv, xmin, x);
                        auto iter = containers.find(key);
                        if (iter == containers.end()) {
                            Container_t container = {{id}};
                            containers.insert({key, container});
                        } else {
                            Container_t &container = iter->second;
                            container.insert({id});
                        }
                    }
                }
            }
        }

        // allocate grid
        return std::unique_ptr<GridSearch>{
            new GridSearch{
                ndim,
                points,
                cells,
                ndiv,
                xmin,
                xmax,
                side_length,
                bounding_boxes,
                containers,
                large_cells,
            }};
    }

    // Checks whether point is in cell or not
    bool in_cell(vector<double> &zeta, vector<double> const &x, CellId id) {
        for (size_t i = 0; i < this->ndim; i++) {
            if (x[i] < this->bounding_boxes[id][i][I_MIN] || x[i] > this->bounding_boxes[id][i][I_MAX]) {
                return false;  // outside the bounding box
            }
        }
        if (this->ndim == 2) {
            triangle_coords(zeta,
                            this->points[this->cells[id][0]],
                            this->points[this->cells[id][1]],
                            this->points[this->cells[id][2]], x);
            return in_triangle(zeta);
        } else {
            tetrahedron_coords(zeta,
                               this->points[this->cells[id][0]],
                               this->points[this->cells[id][1]],
                               this->points[this->cells[id][2]],
                               this->points[this->cells[id][3]], x);
            return in_tetrahedron(zeta);
        }
    }

    // Finds the cell where the given point falls in
    //
    // # Input
    //
    // * `x` -- (ndim) the point coordinates
    //
    // # Output
    //
    // * `status` -- SUCCESS: Indicates that the find function found a cell with success
    //               CELL_NOT_FOUND -- Indicates that the find function could not find a cell
    //               FAIL_POINT_OUTSIDE -- Indicates that the find function failed because
    //                                     the point is outside the grid
    // * `zeta` -- (ndim+1) are the triangle/tetrahedron coordinates
    // * if SUCCESS, returns the index of the cell, otherwise returns -1
    CellId find(size_t &status, vector<double> &zeta, vector<double> const &x) {
        // check
        if (zeta.size() != ndim + 1) {
            throw "the size of zeta must be equal to ndim + 1";
        }

        // check if the point is outside the grid
        for (size_t i = 0; i < ndim; i++) {
            if (x[i] < this->xmin[i] || x[i] > this->xmax[i]) {
                status = FAIL_POINT_OUTSIDE;
                return -1;
            }
        }

        // check if the point is in a large cell
        for (const auto &id : this->large_cells) {
            if (this->in_cell(zeta, x, id)) {
                status = SUCCESS;
                return id;
            }
        }

        // get the container where `x` falls in
        auto key = calc_container_key(this->ndim, this->side_length, this->ndiv, this->xmin, x);
        auto iter = this->containers.find(key);
        if (iter == this->containers.end()) {
            status = CELL_NOT_FOUND;  // there is no container for the key corresponding to x
            return -1;
        }

        // find the cell where the point falls in
        auto container = iter->second;
        for (const auto &id : container) {
            if (this->in_cell(zeta, x, id)) {
                status = SUCCESS;
                return id;
            }
        }

        // not found
        status = CELL_NOT_FOUND;
        return -1;
    }

    // Computes the interpolation after the the cell and zeta are found
    inline double interp(CellId id, vector<double> const &zeta) {
        if (this->ndim == 2) {
            return zeta[0] * this->points[this->cells[id][0]][this->ndim] +
                   zeta[1] * this->points[this->cells[id][1]][this->ndim] +
                   zeta[2] * this->points[this->cells[id][2]][this->ndim];
        } else {
            return zeta[0] * this->points[this->cells[id][0]][this->ndim] +
                   zeta[1] * this->points[this->cells[id][1]][this->ndim] +
                   zeta[2] * this->points[this->cells[id][2]][this->ndim] +
                   zeta[3] * this->points[this->cells[id][3]][this->ndim];
        }
    }

    // Print details about the grid (e.g., for debugging)
    void print_details() {
        cout << "\nGridSearch" << endl;
        cout << "==========" << endl;
        cout << "number of non-empty containers = " << this->containers.size() << endl;
        for (const auto &[key, container] : this->containers) {
            cout << "container # " << key << ": cells = [";
            bool first = true;
            for (const auto &id : container) {
                if (!first) {
                    cout << ", ";
                }
                cout << id;
                first = false;
            }
            cout << "]" << endl;
        }
        cout << "large_cells = [";
        bool first = true;
        for (const auto &id : this->large_cells) {
            if (!first) {
                cout << ", ";
            }
            cout << id;
            first = false;
        }
        cout << "]" << endl;
    }
};
