#pragma once

#include <cmath>
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
bool is_point_inside_triangle(vector<double> const& xa,
                              vector<double> const& xb,
                              vector<double> const& xc,
                              vector<double> const& xp) {
    if ((xp[0] == xa[0] && xp[1] == xa[1]) || (xp[0] == xb[0] && xp[1] == xb[1]) || (xp[0] == xc[0] && xp[1] == xc[1])) {
        return true;
    }
    bool na = cross_is_negative(xb[0] - xa[0], xb[1] - xa[1], xp[0] - xa[0], xp[1] - xa[1]);
    bool nb = cross_is_negative(xc[0] - xb[0], xc[1] - xb[1], xp[0] - xb[0], xp[1] - xb[1]);
    bool nc = cross_is_negative(xa[0] - xc[0], xa[1] - xc[1], xp[0] - xc[0], xp[1] - xc[1]);
    return ((na && nb && nc) || (!na && !nb && !nc));
}
