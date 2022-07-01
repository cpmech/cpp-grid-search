#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

using namespace std;

// Returns true if the sign of the cross product vector is negative
//
// Specifically, returns true if the sign of the component of the out-of-plane vector,
// resulting from the cross product between u and v is negative
inline bool cross_is_negative(double u0, double u1, double v0, double v1) {
    return (u0 * v1 - u1 * v0) < 0.0;
}

// Checks whether a point is inside a triangle or not
//
// ```text
// →   →    →    →   →    →    →   →    →
// a = AB × AP;  b = BC × BP;  c = CA × CP
//
// The point is inside iff:
//     sign(a[2]) == sign(b[2]) == sign(c[2])
//
// A,                      A, -------P   OUTSIDE
// |\',                    | ',     / |
// | \ ',                  |   ',  /    .
// |  \  ',                |     ',      .
// |   \   ',              |     / ',     .
// |    P.   ',            |    /    ',    |
// |   /  '.   ',          |   /       ',   .
// |  /     `-.  ',        |  /          ',  |
// | /  INSIDE `-. ',      | /             ', .
// |/             `-.',    |/                ',|
// C-------------------B   C-------------------B
// ```
//
// **Note:** This function returns true if the point is
// on any boundary or coincides with any vertex
inline bool in_triangle(vector<double> const &xa,
                        vector<double> const &xb,
                        vector<double> const &xc,
                        vector<double> const &xp) {
    if ((xp[0] == xa[0] && xp[1] == xa[1]) || (xp[0] == xb[0] && xp[1] == xb[1]) || (xp[0] == xc[0] && xp[1] == xc[1])) {
        return true;
    }
    bool na = cross_is_negative(xb[0] - xa[0], xb[1] - xa[1], xp[0] - xa[0], xp[1] - xa[1]);
    bool nb = cross_is_negative(xc[0] - xb[0], xc[1] - xb[1], xp[0] - xb[0], xp[1] - xb[1]);
    bool nc = cross_is_negative(xa[0] - xc[0], xa[1] - xc[1], xp[0] - xc[0], xp[1] - xc[1]);
    return ((na && nb && nc) || (!na && !nb && !nc));
}

// GridSearch tolerance for all directions
const double GS_TOLERANCE = 1e-4;

// GridSearch border tolerance to handle imprecision near the borders
const double GS_BORDER_TOL = 1e-2;

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

// Some constants
const size_t NDIM = 2;   // 2D
const size_t NNODE = 3;  // Triangles

// Defines a tool to search the cell where a point is located within a mesh
//
// # Reference
//
// * Durand, Farias, and Pedroso (2015) Computing intersections between
//   non-compatible curves and finite elements, Computational Mechanics;
//   DOI=10.1007/s00466-015-1181-y
struct GridSearch {
    vector<size_t> ndiv;                // (NDIM) number of divisions along each direction
    vector<double> xmin;                // (NDIM) min values
    vector<double> xmax;                // (NDIM) max values
    double side_length;                 // side length of a container
    vector<BboxMinMax> bounding_boxes;  // (NCELL) bounding boxes
    Containers_t containers;            // holds all items

    // Calculates the key of the container where the point should fall in
    static inline ContainerKey calc_container_key(
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
        return ix + iy * ndiv[0];
    }

    // Allocates new instance
    ///
    /// # Input
    ///
    /// * `ncell` -- is the number of cells (e.g., triangle/tetrahedron) in the mesh.
    ///     - All cells are numbered from `0` to `ncell - 1`
    ///     - The index of the cell in a mesh is also called CellId (`cell_id`)
    /// * `get_x` -- is a function of the `cell_id` and the local index of the node/point `m`.
    ///    This function returns the coordinates `x` of the point.
    static std::unique_ptr<GridSearch> make_new(size_t ncell, function<vector<double> const *(size_t, size_t)> get_x) {
        // allocate some variables
        auto xmin = vector<double>(NDIM);
        auto xmax = vector<double>(NDIM);
        auto x_min_max = BboxMinMax(NDIM);
        auto bbox_large = vector<double>(NDIM);
        vector<BboxMinMax> bounding_boxes;
        for (size_t i = 0; i < NDIM; i++) {
            x_min_max[i] = vector<double>(N_MIN_MAX);
            bbox_large[i] = numeric_limits<double>::min();
        }

        // find limits, bounding boxes, and largest cell
        for (size_t cell_id = 0; cell_id < ncell; cell_id++) {
            for (size_t m = 0; m < NNODE; m++) {
                auto x = get_x(cell_id, m);
                if (x->size() != NDIM) {
                    throw "x.size() must be equal to ndim";
                }
                for (size_t i = 0; i < NDIM; i++) {
                    // limits
                    xmin[i] = min(xmin[i], (*x)[i]);
                    xmax[i] = max(xmax[i], (*x)[i]);
                    // bounding box
                    if (m == 0) {
                        for (size_t i = 0; i < NDIM; i++) {
                            x_min_max[i][I_MIN] = (*x)[i];
                            x_min_max[i][I_MAX] = (*x)[i];
                        }
                    } else {
                        x_min_max[i][I_MIN] = min(x_min_max[i][I_MIN], (*x)[i]);
                        x_min_max[i][I_MAX] = max(x_min_max[i][I_MAX], (*x)[i]);
                    }
                }
            }
            // largest cell
            for (size_t i = 0; i < NDIM; i++) {
                bbox_large[i] = max(bbox_large[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            // add to bounding box maps
            bounding_boxes.push_back(x_min_max);
        }

        // make the side_length equal to the largest bounding box dimension
        double side_length = bbox_large[0];
        for (size_t i = 1; i < NDIM; i++) {
            side_length = max(side_length, bbox_large[i]);
        }

        // expand side_length by two times the tolerance (left and right)
        side_length += 2.0 * GS_TOLERANCE;

        // expand borders
        for (size_t i = 0; i < NDIM; i++) {
            xmin[i] -= GS_BORDER_TOL;
            xmax[i] += GS_BORDER_TOL;
        }

        // number of divisions
        vector<size_t> ndiv(NDIM);
        for (size_t i = 0; i < NDIM; i++) {
            ndiv[i] = ceil((xmax[i] - xmin[i]) / side_length);
        }

        // update xmax after deciding on the side_length and number of divisions
        for (size_t i = 0; i < NDIM; i++) {
            xmax[i] = xmin[i] + side_length * ((double)ndiv[i]);
        }

        // insert cells
        Containers_t containers;
        vector<double> x(NDIM);
        for (size_t cell_id = 0; cell_id < ncell; cell_id++) {
            auto x_min_max = bounding_boxes[cell_id];
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    x[0] = x_min_max[0][r];
                    x[1] = x_min_max[1][s];
                    auto key = calc_container_key(side_length, ndiv, xmin, x);
                    auto iter = containers.find(key);
                    if (iter == containers.end()) {
                        Container_t container = {{cell_id}};
                        containers.insert({key, container});
                    } else {
                        Container_t &container = iter->second;
                        container.insert({cell_id});
                    }
                }
            }
        }

        // allocate grid
        return std::unique_ptr<GridSearch>{
            new GridSearch{
                ndiv,
                xmin,
                xmax,
                side_length,
                bounding_boxes,
                containers,
            }};
    }

    // Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
    //
    // * Returns the CellId or -1 if no cell contains the point
    // * `is_in_cell` is a function of cell_id and x that tells whether the point si in the cell or not
    int find_cell(vector<double> &x, function<bool(size_t, vector<double> const *)> is_in_cell) {
        // check if the point is out-of-bounds
        for (size_t i = 0; i < NDIM; i++) {
            if (x[i] < this->xmin[i] || x[i] > this->xmax[i]) {
                throw "given point coordinates are outsize the grid";
            }
        }

        // get the container where `x` falls in
        auto key = calc_container_key(this->side_length, this->ndiv, this->xmin, x);
        auto iter = this->containers.find(key);
        if (iter == this->containers.end()) {
            return -1;  // there is no container for the key corresponding to x
        }

        // find the cell where the point falls in
        auto container = iter->second;
        for (const auto &cell_id : container) {
            auto x_min_max = this->bounding_boxes[cell_id];
            for (size_t i = 0; i < NDIM; i++) {
                if (x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX]) {
                    continue;  // outside the bounding box
                }
            }
            if ((is_in_cell)(cell_id, &x)) {
                return cell_id;
            }
        }

        // not found
        return -1;
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
    }
};
