#pragma once

#include <map>
#include <memory>
#include <set>
#include <vector>

using namespace std;

// Specifies the key of containers (or bins in the grid)
typedef size_t ContainerKey;

// Specifies the identification number of items
typedef size_t ItemID;

// Defines the container type
typedef set<ItemID> Container_t;

// Defines the containers type
typedef map<ContainerKey, Container_t> Containers_t;

// Holds the options for GridSearch
//
// * `xmin` -- (ndim) minimum coordinates to define the boundary
// * `xmax` -- (ndim) maximum coordinates to define the boundary
// * `expansion` -- expansion factor for the bounding box.
//    This constant serves to accommodate eventual imprecision near the boundaries.
//    The value is a multiplier for `delta = max - min`
//    For example, **(recommended) expansion = 0.01** means 1% of delta
// * `ndiv` -- number of division along each direction
// * `tols` -- tolerance for each direction; e.g. [1e-4, 1e-4]
struct GridSearchOptions {
    double expansion;
    vector<double> xmin;
    vector<double> xmax;
    vector<double> tols;
    vector<size_t> ndiv;

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
    size_t ndim;                 // space dimension
    double expansion;            // expansion factor; e.g. 0.01
    vector<double> xmin;         // (ndim) min values
    vector<double> xmax;         // (ndim) max values
    vector<double> tols;         // (ndim) tolerances
    vector<size_t> ndiv;         // (ndim) number of divisions along each direction
    vector<double> delta;        // (ndim) difference between max and min
    vector<double> side_length;  // (ndim) side lengths of each container
    vector<size_t> coefficient;  // (3) coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)
    Containers_t containers;     // holds all items

    // Allocates new instance
    static std::unique_ptr<GridSearch> make_new(const std::unique_ptr<GridSearchOptions> &options);

    // Inserts cells (e.g., triangles or tetrahedra) to the correct containers in the grid
    void insert_cells(vector<vector<vector<double>>> &mesh);

    // Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
    //
    // # Input
    //
    // * `x` -- coordinates (ndim)
    //
    // # Output
    //
    // * `id` -- if found, returns the ID of the cell,
    //           otherwise returns -1 (not found)
    int find_cell(vector<double> &x, vector<vector<vector<double>>> &mesh);

    // Print details about the grid (e.g., for debugging)
    void print_details();

   private:
    // Calculates the key of the container where the point x should be located
    //
    // # Output
    //
    // * returns the key of the container or -1 if the point is out-of-range (i.e., outside [xmin,xmax])
    int _get_container_key(vector<double> &x);

    // Updates a container or inserts a point into an existing container
    void _update_or_insert_point(ContainerKey key, ItemID id, vector<double> &x);
};