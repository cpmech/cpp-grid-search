#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

extern "C" {
#include "triangle/constants.h"
#include "triangle/interface_triangle.h"
}

using namespace std;

// Returns a new Delaunay triangulation from a "cloud" of points
//
// Input:
//
// * `cloud` -- is a list of coordinates such as `[[x0,y0,...], [x1,y1,..], ...]`
// * `verbose` -- prints Triangle's detailed output on the console
//
// Returns the list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
vector<vector<size_t>> delaunay_2d(vector<vector<double>> const& cloud, bool verbose) {
    // constants
    size_t npoint = cloud.size();
    if (npoint < 3) {
        throw "number of points must be >= 3";
    }
    size_t ncol = cloud[0].size();
    if (ncol < 2) {
        throw "coordinates.ncol must be >= 2";
    }

    // allocate triangle data
    auto ext_triangle = new_triangle(npoint, 0, 0, 0);
    if (ext_triangle == NULL) {
        throw "cannot allocate triangle data";
    }

    // set points
    for (size_t index = 0; index < npoint; index++) {
        int status = set_point(ext_triangle, index, cloud[index][0], cloud[index][1]);
        if (status != TRITET_SUCCESS) {
            if (status == TRITET_ERROR_NULL_DATA) {
                throw "INTERNAL ERROR: found NULL data";
            }
            if (status == TRITET_ERROR_NULL_POINT_LIST) {
                throw "INTERNAL ERROR: found NULL point list";
            }
            if (status == TRITET_ERROR_INVALID_POINT_INDEX) {
                throw "index of point is out of bounds";
            }
            throw "INTERNAL ERROR: some error occurred";
        }
    }

    // perform Delaunay triangulation
    int status = run_delaunay(ext_triangle, verbose);
    if (status != TRITET_SUCCESS) {
        if (status == TRITET_ERROR_NULL_DATA) {
            throw "INTERNAL ERROR: found NULL data";
        }
        if (status == TRITET_ERROR_NULL_POINT_LIST) {
            throw "INTERNAL ERROR: found NULL point list";
        }
        throw "INTERNAL ERROR: some error occurred";
    }

    // extract points from Triangle' data structures
    size_t n_mesh_point = get_npoint(ext_triangle);
    auto coordinates = vector<vector<double>>(n_mesh_point);
    for (size_t point_index = 0; point_index < n_mesh_point; point_index++) {
        coordinates[point_index] = vector<double>(3);  // x,y,T (in case we want to use the third value for interpolations)
        coordinates[point_index][0] = get_point(ext_triangle, point_index, 0);
        coordinates[point_index][1] = get_point(ext_triangle, point_index, 1);
        coordinates[point_index][2] = 0.0;  // zero "temperature"
    }

    // extract triangles from Triangle' data structures
    size_t ntriangle = get_ntriangle(ext_triangle);
    size_t ncorner = get_ncorner(ext_triangle);  // this should be 3 because we didn't require "o2" triangles
    auto triangles = vector<vector<size_t>>(ntriangle);
    for (size_t tri_index = 0; tri_index < ntriangle; tri_index++) {
        triangles[tri_index] = vector<size_t>(ncorner);
        for (size_t corner = 0; corner < ncorner; corner++) {
            triangles[tri_index][corner] = get_triangle_corner(ext_triangle, tri_index, corner);
        }
    }

    // make sure we clean Triangle data
    drop_triangle(ext_triangle);

    // results
    return triangles;
}
