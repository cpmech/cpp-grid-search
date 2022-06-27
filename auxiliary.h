#pragma once

#include <cmath>
#include <vector>

using namespace std;

// Computes the point-to-point distance
//
// # Input
//
// * `a` -- first point
// * `b` -- second point
//
// # Output
//
// Returns the unsigned distance between the two points.
inline double point_point_distance(vector<double> const &a, vector<double> const &b) {
    size_t ndim = a.size();
    if (ndim < 2 || ndim > 3) {
        throw "a.size() == ndim must be 2 or 3";
    }
    if (b.size() != ndim) {
        throw "b.size() must equal a.size() == ndim";
    }
    if (ndim == 2) {
        return sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));
    } else {
        return sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]));
    }
}

inline double tri_sign(vector<double> const &a, vector<double> const &b, vector<double> const &c) {
    // From here: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
    return (a[0] - c[0]) * (b[1] - c[1]) - (b[0] - c[0]) * (a[1] - c[1]);
}

inline bool is_point_in_triangle(vector<double> const &p, vector<vector<double>> const &points) {
    // From here: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle

    double d1 = tri_sign(p, points[0], points[1]);
    double d2 = tri_sign(p, points[1], points[2]);
    double d3 = tri_sign(p, points[2], points[0]);

    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

inline bool is_point_in_tetrahedron(vector<double> const &p, vector<vector<double>> const &points) {
    // TODO
    return false;
}

inline void compute_bounding_box(double expansion, vector<double> &xmin, vector<double> &xmax, vector<vector<double>> const &cell) {
    size_t ndim = xmin.size();
    if (ndim < 2 || ndim > 3) {
        throw "xmin.size() == ndim must be 2 or 3";
    }
    if (xmax.size() != ndim) {
        throw "xmax.size() must equal xmin.size() == ndim";
    }
    size_t nnode = cell.size();
    if (nnode < 3 || nnode > 4) {
        throw "cell.size() == nnode must be 3 (triangle) or 4 (tetrahedron)";
    }
    for (size_t i = 0; i < ndim; i++) {
        xmin[i] = cell[0][i];
        xmax[i] = cell[0][i];
        for (size_t m = 1; m < nnode; m++) {
            xmin[i] = min(xmin[i], cell[m][i]);
            xmax[i] = max(xmax[i], cell[m][i]);
        }
    }
    for (size_t i = 0; i < ndim; i++) {
        double delta = xmax[i] - xmin[i];
        xmin[i] -= expansion * delta;
        xmax[i] += expansion * delta;
    }
}
