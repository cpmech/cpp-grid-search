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