#pragma once

#include <vector>

#include "tetgen/interface_tetgen.h"
#include "triangle/constants.h"

using namespace std;

// Returns a new Delaunay tetrahedralization from a "cloud" of points
//
// Input:
//
// * `cloud` -- is a list of coordinates such as `[[x0,y0,z0,...], [x1,y1,z1,..], ...]`
// * `verbose` -- prints TetGen's detailed output on the console
//
// Returns the list of connectivity (topology) such as `[[0,2,1,3], [2,1,0,4]]`
vector<vector<size_t>> delaunay_3d(vector<vector<double>> const& cloud, bool verbose) {
    // constants
    size_t npoint = cloud.size();
    if (npoint < 4) {
        throw "number of points must be >= 4";
    }
    size_t ncol = cloud[0].size();
    if (ncol < 3) {
        throw "coordinates.ncol must be >= 3";
    }

    // allocate Tetgen data
    auto ext_tetgen = new_tetgen(npoint, 0, 0, 0, 0);
    if (ext_tetgen == NULL) {
        throw "cannot allocate tetgen data";
    }

    // set points
    for (size_t index = 0; index < npoint; index++) {
        int status = tet_set_point(ext_tetgen, index, cloud[index][0], cloud[index][1], cloud[index][2]);
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
    int status = tet_run_delaunay(ext_tetgen, verbose);
    if (status != TRITET_SUCCESS) {
        if (status == TRITET_ERROR_NULL_DATA) {
            throw "INTERNAL ERROR: found NULL data";
        }
        if (status == TRITET_ERROR_NULL_POINT_LIST) {
            throw "INTERNAL ERROR: found NULL point list";
        }
        throw "INTERNAL ERROR: some error occurred";
    }

    // extract points from Tetgen' data structures
    size_t n_mesh_point = tet_get_npoint(ext_tetgen);
    auto coordinates = vector<vector<double>>(n_mesh_point);
    for (size_t point_index = 0; point_index < n_mesh_point; point_index++) {
        coordinates[point_index] = vector<double>(4);  // x,y,z,T (in case we want to use the fourth value for interpolations)
        coordinates[point_index][0] = tet_get_point(ext_tetgen, point_index, 0);
        coordinates[point_index][1] = tet_get_point(ext_tetgen, point_index, 1);
        coordinates[point_index][2] = tet_get_point(ext_tetgen, point_index, 2);
        coordinates[point_index][3] = 0.0;  // zero "temperature"
    }

    // extract tetrahedra from Tetgen' data structures
    size_t ntet = tet_get_ntetrahedron(ext_tetgen);
    size_t ncorner = tet_get_ncorner(ext_tetgen);  // this should be 4 because we didn't require "o2" tetrahedra
    auto tetrahedra = vector<vector<size_t>>(ntet);
    for (size_t tet_index = 0; tet_index < ntet; tet_index++) {
        tetrahedra[tet_index] = vector<size_t>(ncorner);
        for (size_t corner = 0; corner < ncorner; corner++) {
            tetrahedra[tet_index][corner] = tet_get_tetrahedron_corner(ext_tetgen, tet_index, corner);
        }
    }

    // make sure we clean Tetgen data
    drop_tetgen(ext_tetgen);

    // results
    return tetrahedra;
}
