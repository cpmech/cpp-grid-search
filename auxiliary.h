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
inline double point_point_distance(vector<double> &a, vector<double> &b) {
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

inline double tri_sign(vector<double> &a, vector<double> &b, vector<double> &c) {
    // From here: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
    return (a[0] - c[0]) * (b[1] - c[1]) - (b[0] - c[0]) * (a[1] - c[1]);
}

inline bool is_point_in_triangle(vector<double> &p, vector<vector<double>> points) {
    // From here: https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle

    double d1 = tri_sign(p, points[0], points[1]);
    double d2 = tri_sign(p, points[1], points[2]);
    double d3 = tri_sign(p, points[2], points[0]);

    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}