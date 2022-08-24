#pragma once

#include <vector>

using namespace std;

// Calculates the triangle (three) coordinates
//
// # Output
//
// * `zeta` -- the triangle coordinates (len = 3)
//
// # Input
//
// * `xa,xb,xc` -- corners (each has len = 2)
// * `xp` -- the point to calculate the coordinates (len = 2)
//
// # Exceptions
//
// This function will throw exception if the array sizes are incorrect
inline void triangle_coords(vector<double> &zeta,
                            vector<double> const &xa,
                            vector<double> const &xb,
                            vector<double> const &xc,
                            vector<double> const &xp) {
    double a2 = xa[0] * (xb[1] - xc[1]) + xb[0] * (xc[1] - xa[1]) + xc[0] * (xa[1] - xb[1]);
    zeta[0] = (xp[0] * (xb[1] - xc[1]) + xb[0] * (xc[1] - xp[1]) + xc[0] * (xp[1] - xb[1])) / a2;
    zeta[1] = (xa[0] * (xp[1] - xc[1]) + xp[0] * (xc[1] - xa[1]) + xc[0] * (xa[1] - xp[1])) / a2;
    zeta[2] = (xa[0] * (xb[1] - xp[1]) + xb[0] * (xp[1] - xa[1]) + xp[0] * (xa[1] - xb[1])) / a2;
}

// Indicates if a point is inside a triangle by looking at its triangle coordinates (zeta)
//
// Note: the point is inside (or on an edge) if all zeta are positive (or zero)
//
// # Input
//
// * `zeta` -- the triangle coordinates (len = 3)
//
// # Exceptions
//
// This function will throw exception if the array sizes are incorrect
inline bool in_triangle(vector<double> const &zeta) {
    if (zeta[0] < 0.0 || zeta[1] < 0.0 || zeta[2] < 0.0) {
        return false;
    } else {
        return true;
    }
}

// Calculates the tetrahedron (four) coordinates
//
// # Output
//
// * `zeta` -- the tetrahedron coordinates (len = 4)
//
// # Input
//
// * `xa,xb,xc,xd` -- corners (each has len = 3)
// * `xp` -- the point to calculate the coordinates (len = 3)
//
// # Exceptions
//
// This function will throw exception if the array sizes are incorrect
inline void tetrahedron_coords(vector<double> &zeta,
                               vector<double> const &xa,
                               vector<double> const &xb,
                               vector<double> const &xc,
                               vector<double> const &xd,
                               vector<double> const &xp) {
    double x1 = xa[0];
    double y1 = xa[1];
    double z1 = xa[2];
    double x2 = xb[0];
    double y2 = xb[1];
    double z2 = xb[2];
    double x3 = xc[0];
    double y3 = xc[1];
    double z3 = xc[2];
    double x4 = xd[0];
    double y4 = xd[1];
    double z4 = xd[2];
    double x = xp[0];
    double y = xp[1];
    double z = xp[2];
    double v6 = (x2 - x1) * ((y2 - y3) * (z3 - z4) - (y3 - y4) * (z2 - z3)) + (x3 - x2) * ((y3 - y4) * (z1 - z2) - (y1 - y2) * (z3 - z4)) + (x4 - x3) * ((y1 - y2) * (z2 - z3) - (y2 - y3) * (z1 - z2));
    zeta[0] = ((x2 * (y3 * z4 - y4 * z3) + x3 * (y4 * z2 - y2 * z4) + x4 * (y2 * z3 - y3 * z2)) + x * ((y4 - y2) * (z3 - z2) - (y3 - y2) * (z4 - z2)) + y * ((x3 - x2) * (z4 - z2) - (x4 - x2) * (z3 - z2)) + z * ((x4 - x2) * (y3 - y2) - (x3 - x2) * (y4 - y2))) / v6;
    zeta[1] = ((x1 * (y4 * z3 - y3 * z4) + x3 * (y1 * z4 - y4 * z1) + x4 * (y3 * z1 - y1 * z3)) + x * ((y3 - y1) * (z4 - z3) - (y3 - y4) * (z1 - z3)) + y * ((x4 - x3) * (z3 - z1) - (x1 - x3) * (z3 - z4)) + z * ((x3 - x1) * (y4 - y3) - (x3 - x4) * (y1 - y3))) / v6;
    zeta[2] = ((x1 * (y2 * z4 - y4 * z2) + x2 * (y4 * z1 - y1 * z4) + x4 * (y1 * z2 - y2 * z1)) + x * ((y2 - y4) * (z1 - z4) - (y1 - y4) * (z2 - z4)) + y * ((x1 - x4) * (z2 - z4) - (x2 - x4) * (z1 - z4)) + z * ((x2 - x4) * (y1 - y4) - (x1 - x4) * (y2 - y4))) / v6;
    zeta[3] = ((x1 * (y3 * z2 - y2 * z3) + x2 * (y1 * z3 - y3 * z1) + x3 * (y2 * z1 - y1 * z2)) + x * ((y1 - y3) * (z2 - z1) - (y1 - y2) * (z3 - z1)) + y * ((x2 - x1) * (z1 - z3) - (x3 - x1) * (z1 - z2)) + z * ((x1 - x3) * (y2 - y1) - (x1 - x2) * (y3 - y1))) / v6;
}

// Indicates if a point is inside a tetrahedron by looking at its tetrahedron coordinates (zeta)
//
// Note: the point is inside (or on a face) if all zeta are positive (or zero)
//
// # Input
//
// * `zeta` -- the tetrahedron coordinates (len = 4)
//
// # Exceptions
//
// This function will throw exception if the array sizes are incorrect
inline bool in_tetrahedron(vector<double> const &zeta) {
    if (zeta[0] < 0.0 || zeta[1] < 0.0 || zeta[2] < 0.0 || zeta[3] < 0.0) {
        return false;
    } else {
        return true;
    }
}
