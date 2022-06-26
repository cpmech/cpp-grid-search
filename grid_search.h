#pragma once

#include <map>
#include <memory>
#include <vector>

using namespace std;

struct Item {
    size_t id;         // identification number
    vector<double> x;  // (ndim) coordinates
};

// Holds the options for GridSearch
//
// * `min` -- (ndim) minimum coordinates to define the boundary
// * `max` -- (ndim) maximum coordinates to define the boundary
// * `expansion` -- expansion factor for the bounding box.
//    This constant serves to accommodate eventual imprecision near the boundaries.
//    The value is a multiplier for `delta = max - min`
//    For example, **(recommended) expansion = 0.01** means 1% of delta
// * `ndiv` -- number of division along each direction
// * `tol` -- tolerance; e.g. 1e-4
struct GridSearchOptions {
    vector<double> min;
    vector<double> max;
    double expansion;
    vector<size_t> ndiv;
    double tol;

    // Allocates new instance
    static std::unique_ptr<GridSearchOptions> make_new();
};

// Implements the GridSearch tool
struct GridSearch {
    // constants
    size_t ndim;           // space dimension
    vector<size_t> ndiv;   // (ndim) number of divisions along each direction
    vector<double> min;    // (ndim) min values
    vector<double> max;    // (ndim) max values
    vector<double> delta;  // (ndim) difference between max and min
    vector<double> size;   // (ndim) side lengths of each container
    vector<size_t> cf;     // (3) coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)

    // square/cubic halo: bounding box corners around point, including the point
    vector<double> tol;           // tolerances to compare coordinates and define the halo
    vector<vector<double>> halo;  // (ncorner) 4 in 2D or 8 in 3D (each contains ndim coords)
    size_t ncorner;               // number of halo corners 4 in 2D or 8 in 3D

    // holds non-empty containers. maps container.index to container.data
    // a point may be located in more than one container (e.g., when at internal boundaries)
    map<size_t, map<size_t, Item>> containers;

    // constants
    double radius;      // radius of the circumscribed circle of containers
    double radius_tol;  // radius of the bounding box defined by the tolerances

    // Allocates new instance
    static std::unique_ptr<GridSearch> make_new(const std::unique_ptr<GridSearchOptions> &options);
};