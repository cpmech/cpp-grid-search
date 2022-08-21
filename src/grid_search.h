#pragma once

#include <cmath>
#include <iostream>
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

// Interpolates a value at a given point over a triangle
//
// See file data/derivations/interpolation-over-triangle.pdf
//
// # Input
//
// `xa,xb,xc` -- are the coordinates of the triangle vertices (each: len = 2 = ndim)
// `temp` -- are the known values at the triangle vertices; e.g., the temperatures (len = 3)
// `xp` -- are the coordinates of the point of interest (len = ndim)
//
// # Output
//
// Returns the "temperature" `T @ x`
//
// # Panics
//
// This function will panic if the array sizes are incorrect
inline double triangle_interpolation(
    vector<double> const &xa,
    vector<double> const &xb,
    vector<double> const &xc,
    vector<double> const &temp,
    vector<double> const &xp) {
    double aa2 = xa[0] * (xb[1] - xc[1]) + xb[0] * (xc[1] - xa[1]) + xc[0] * (xa[1] - xb[1]);
    double zeta_0 = (xp[0] * (xb[1] - xc[1]) + xb[0] * (xc[1] - xp[1]) + xc[0] * (xp[1] - xb[1])) / aa2;
    double zeta_1 = (xa[0] * (xp[1] - xc[1]) + xp[0] * (xc[1] - xa[1]) + xc[0] * (xa[1] - xp[1])) / aa2;
    double zeta_2 = (xa[0] * (xb[1] - xp[1]) + xb[0] * (xp[1] - xa[1]) + xp[0] * (xa[1] - xb[1])) / aa2;
    return zeta_0 * temp[0] + zeta_1 * temp[1] + zeta_2 * temp[2];
}

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
/// The other constant [GS_LARGE_CELLS_MAX_RATIO] also influences the algorithm by enabling or disabling this strategy.
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

// Specifies the identification number of triangles (must be sequential from 0 to ntriangle - 1)
typedef size_t TriangleId;

// Defines the container type
typedef set<TriangleId> Container_t;

// Defines the containers type
typedef map<ContainerKey, Container_t> Containers_t;

// Defines the bounding box of a triangle
const size_t N_MIN_MAX = 2;  // 2 means {min,max}
const size_t I_MIN = 0;
const size_t I_MAX = 1;
typedef vector<vector<double>> BboxMinMax;  // [ndim][N_MIN_MAX]

// Some constants
const size_t NDIM = 2;   // 2D
const size_t NNODE = 3;  // Triangles

// Defines a tool to search the triangle where a point is located within a mesh
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
    vector<BboxMinMax> bounding_boxes;  // (NTRIANGLE) bounding boxes
    Containers_t containers;            // holds all items
    Container_t large_triangles;        // holds the ids of large triangles

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
    //
    // # Input
    //
    // * `coordinates` -- is a list of coordinates such as `[[x0,y0,T0], [x1,y1,T1], [x2,y2,T2], [x3,y3,T3]]`
    //   where `T[i]` are the values (e.g., temperatures) at that coordinate and will be used for interpolation
    // * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    static std::unique_ptr<GridSearch> make_new(vector<vector<double>> const &coordinates,
                                                vector<vector<size_t>> const &triangles) {
        // constants
        size_t npoint = coordinates.size();
        if (npoint < 3) {
            throw "number of points must be >= 3";
        }
        size_t ncol = coordinates[0].size();
        if (ncol < 2) {
            throw "coordinates.ncol must be >= 2";
        }
        size_t ntriangle = triangles.size();
        if (ntriangle < 1) {
            throw "number of triangles must be >= 1";
        }
        size_t nnode = triangles[0].size();
        if (nnode != 3) {
            throw "number of triangle nodes (nnode) must be = 3";
        }

        // allocate some variables
        auto x = vector<double>(NDIM);
        auto xmin = vector<double>(NDIM);
        auto xmax = vector<double>(NDIM);
        auto x_min_max = BboxMinMax(NDIM);
        auto bbox_largest = vector<double>(NDIM);
        vector<BboxMinMax> bounding_boxes;
        for (size_t i = 0; i < NDIM; i++) {
            xmin[i] = numeric_limits<double>::max();
            xmax[i] = numeric_limits<double>::min();
            x_min_max[i] = vector<double>(N_MIN_MAX);
            bbox_largest[i] = numeric_limits<double>::min();
        }

        // find limits, bounding boxes, and largest triangle
        for (size_t id = 0; id < ntriangle; id++) {
            for (size_t m = 0; m < NNODE; m++) {
                size_t p = triangles[id][m];
                x[0] = coordinates[p][0];
                x[1] = coordinates[p][1];
                for (size_t i = 0; i < NDIM; i++) {
                    // limits
                    xmin[i] = min(xmin[i], x[i]);
                    xmax[i] = max(xmax[i], x[i]);
                    // bounding box
                    if (m == 0) {
                        for (size_t i = 0; i < NDIM; i++) {
                            x_min_max[i][I_MIN] = x[i];
                            x_min_max[i][I_MAX] = x[i];
                        }
                    } else {
                        x_min_max[i][I_MIN] = min(x_min_max[i][I_MIN], x[i]);
                        x_min_max[i][I_MAX] = max(x_min_max[i][I_MAX], x[i]);
                    }
                }
            }
            // largest triangle
            for (size_t i = 0; i < NDIM; i++) {
                bbox_largest[i] = max(bbox_largest[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            // add to bounding box maps
            bounding_boxes.push_back(x_min_max);
        }

        // compute largest length of largest bounding box
        double bbox_largest_length = bbox_largest[0];
        for (size_t i = 1; i < NDIM; i++) {
            bbox_largest_length = max(bbox_largest_length, bbox_largest[i]);
        }

        // handle large triangles
        double bbox_large_length = 0.0;  // largest length among the not very large triangles
        Container_t large_triangles;
        for (size_t id = 0; id < ntriangle; id++) {
            auto x_min_max = bounding_boxes[id];
            double tri_largest_length = x_min_max[0][I_MAX] - x_min_max[0][I_MIN];
            for (size_t i = 0; i < NDIM; i++) {
                tri_largest_length = max(tri_largest_length, x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            if (tri_largest_length >= GS_LARGE_CELLS_CUTOFF * bbox_largest_length) {
                large_triangles.insert(id);
            } else {
                for (size_t i = 0; i < NDIM; i++) {
                    bbox_large_length = max(bbox_large_length, x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
                }
            }
        }

        // make the side_length equal to the largest bounding box dimension
        double pct = ((double)(large_triangles.size())) / ((double)(ntriangle));
        double side_length;
        if (pct <= GS_LARGE_CELLS_MAX_COUNT_PCT) {
            side_length = bbox_large_length;  // use the largest length among the not very large triangles
        } else {
            large_triangles.clear();            // abandon the large triangles strategy
            side_length = bbox_largest_length;  // use the largest length among them all
        };

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

        // insert triangles
        Containers_t containers;
        for (size_t id = 0; id < ntriangle; id++) {
            auto iter = large_triangles.find(id);
            if (iter != large_triangles.end()) {
                continue;  // skip large triangles
            }
            auto x_min_max = bounding_boxes[id];
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    x[0] = x_min_max[0][r];
                    x[1] = x_min_max[1][s];
                    auto key = calc_container_key(side_length, ndiv, xmin, x);
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

        // allocate grid
        return std::unique_ptr<GridSearch>{
            new GridSearch{
                ndiv,
                xmin,
                xmax,
                side_length,
                bounding_boxes,
                containers,
                large_triangles,
            }};
    }

    // Find the triangle where the given coordinate falls in
    //
    // # Input
    //
    // * `coordinates` -- is a list of coordinates such as `[[x0,y0], [x1,y1], [x2,y2], [x3,y3]]`
    // * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    //
    // # Output
    //
    // Returns the index of the triangle in `triangles` or -1 if no triangle contains the point
    //
    // # Warning (Exceptions)
    //
    // The pair `coordinates` and `triangles` must be the same as the ones used in the `new` function,
    // otherwise **panics** may occur or, even worse, you may get **incorrect results**.
    int find_triangle(vector<double> &x,
                      vector<vector<double>> const &coordinates,
                      vector<vector<size_t>> const &triangles) {
        // check if the point is out-of-bounds
        for (size_t i = 0; i < NDIM; i++) {
            if (x[i] < this->xmin[i] || x[i] > this->xmax[i]) {
                throw "given point coordinates are outside the grid";
            }
        }

        // auxiliary
        vector<double> xa(NDIM);
        vector<double> xb(NDIM);
        vector<double> xc(NDIM);

        // check if the point is in a large triangle
        for (const auto &id : this->large_triangles) {
            size_t a = triangles[id][0];
            size_t b = triangles[id][1];
            size_t c = triangles[id][2];
            auto x_min_max = this->bounding_boxes[id];
            bool outside_bbox = false;
            for (size_t i = 0; i < NDIM; i++) {
                if (x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX]) {
                    outside_bbox = true;  // outside the bounding box
                    break;
                }
                xa[i] = coordinates[a][i];
                xb[i] = coordinates[b][i];
                xc[i] = coordinates[c][i];
            }
            if (outside_bbox) {
                continue;  // don't even bother calling is_point_inside_triangle
            }
            if (in_triangle(xa, xb, xc, x)) {
                return id;
            }
        }

        // get the container where `x` falls in
        auto key = calc_container_key(this->side_length, this->ndiv, this->xmin, x);
        auto iter = this->containers.find(key);
        if (iter == this->containers.end()) {
            return -1;  // there is no container for the key corresponding to x
        }

        // find the triangle where the point falls in
        auto container = iter->second;
        for (const auto &id : container) {
            size_t a = triangles[id][0];
            size_t b = triangles[id][1];
            size_t c = triangles[id][2];
            auto x_min_max = this->bounding_boxes[id];
            bool outside_bbox = false;
            for (size_t i = 0; i < NDIM; i++) {
                if (x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX]) {
                    outside_bbox = true;  // outside the bounding box
                    break;
                }
                xa[i] = coordinates[a][i];
                xb[i] = coordinates[b][i];
                xc[i] = coordinates[c][i];
            }
            if (outside_bbox) {
                continue;  // don't even bother calling is_point_inside_triangle
            }
            if (in_triangle(xa, xb, xc, x)) {
                return id;
            }
        }

        // not found
        return -1;
    }

    // Find the triangle where the given coordinate falls in and interpolate coordinates
    //
    // # Input
    //
    // * `coordinates` -- is a list of coordinates such as `[[x0,y0], [x1,y1], [x2,y2], [x3,y3]]`
    // * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    //
    // # Output
    //
    // Returns the value (e.g., temperature) at the target point (xp) inside the triangle.
    // Returns NaN if no triangle contains the point.
    //
    // # Warning (Exceptions)
    //
    // The pair `coordinates` and `triangles` must be the same as the ones used in the `new` function,
    // otherwise **panics** may occur or, even worse, you may get **incorrect results**.
    double find_triangle_and_interpolate(vector<double> &x,
                                         vector<vector<double>> const &coordinates,
                                         vector<vector<size_t>> const &triangles) {
        // check if the point is out-of-bounds
        for (size_t i = 0; i < NDIM; i++) {
            if (x[i] < this->xmin[i] || x[i] > this->xmax[i]) {
                throw "given point coordinates are outside the grid";
            }
        }

        // check if the temperature is present in the coordinates list
        if (coordinates[0].size() < 3) {
            throw "coordinates must contain a third column with the temperature values";
        }

        // auxiliary
        vector<double> xa(NDIM);
        vector<double> xb(NDIM);
        vector<double> xc(NDIM);
        vector<double> temp(3);  // 3 nodes

        // check if the point is in a large triangle
        for (const auto &id : this->large_triangles) {
            size_t a = triangles[id][0];
            size_t b = triangles[id][1];
            size_t c = triangles[id][2];
            auto x_min_max = this->bounding_boxes[id];
            bool outside_bbox = false;
            for (size_t i = 0; i < NDIM; i++) {
                if (x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX]) {
                    outside_bbox = true;  // outside the bounding box
                    break;
                }
                xa[i] = coordinates[a][i];
                xb[i] = coordinates[b][i];
                xc[i] = coordinates[c][i];
            }
            if (outside_bbox) {
                continue;
            }
            if (in_triangle(xa, xb, xc, x)) {
                temp[0] = coordinates[a][2];
                temp[1] = coordinates[b][2];
                temp[2] = coordinates[c][2];
                return triangle_interpolation(xa, xb, xc, temp, x);
            }
        }

        // get the container where `x` falls in
        auto key = calc_container_key(this->side_length, this->ndiv, this->xmin, x);
        auto iter = this->containers.find(key);
        if (iter == this->containers.end()) {
            return std::numeric_limits<double>::quiet_NaN();  // there is no container for the key corresponding to x
        }

        // find the triangle where the point falls in
        auto container = iter->second;
        for (const auto &id : container) {
            size_t a = triangles[id][0];
            size_t b = triangles[id][1];
            size_t c = triangles[id][2];
            auto x_min_max = this->bounding_boxes[id];
            bool outside_bbox = false;
            for (size_t i = 0; i < NDIM; i++) {
                if (x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX]) {
                    outside_bbox = true;  // outside the bounding box
                    break;
                }
                xa[i] = coordinates[a][i];
                xb[i] = coordinates[b][i];
                xc[i] = coordinates[c][i];
            }
            if (outside_bbox) {
                continue;
            }
            if (in_triangle(xa, xb, xc, x)) {
                temp[0] = coordinates[a][2];
                temp[1] = coordinates[b][2];
                temp[2] = coordinates[c][2];
                return triangle_interpolation(xa, xb, xc, temp, x);
            }
        }

        // not found
        return std::numeric_limits<double>::quiet_NaN();
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
        cout << "large_triangles = [";
        bool first = true;
        for (const auto &id : this->large_triangles) {
            if (!first) {
                cout << ", ";
            }
            cout << id;
            first = false;
        }
        cout << "]" << endl;
    }
};
