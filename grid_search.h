#pragma once

#include <map>
#include <memory>
#include <vector>

using namespace std;

// Specifies the index of containers (or bins in the grid)
typedef size_t Index;

// Specifies the identification number of items
typedef size_t ID;

// Holds the id and coordinates of an item
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

// Implements a grid for fast searching entries by coordinates
//
// # Reference
//
// * Durand, Farias, and Pedroso (2015) Computing intersections between
//   non-compatible curves and finite elements, Computational Mechanics;
//   DOI=10.1007/s00466-015-1181-y
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
    map<Index, map<ID, Item>> containers;

    // constants
    double radius;      // radius of the circumscribed circle of containers
    double radius_tol;  // radius of the bounding box defined by the tolerances

    // Allocates new instance
    static std::unique_ptr<GridSearch> make_new(const std::unique_ptr<GridSearchOptions> &options);

    // Inserts a new item to the right container in the grid
    void insert(size_t id, vector<double> &x);

    // Calculates the container index where the point x should be located
    int container_index(vector<double> &x);

    // Updates container or inserts point in an existing container
    void update_or_insert(Index index, ID id, vector<double> &x);

    // Sets square/cubic halo around point
    void set_halo(vector<double> &x);

    // Print details about the grid (e.g., for debugging)
    void print_details();
};