#include "grid_search.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace std;

// Allocates new instance
std::unique_ptr<GridSearchOptions> GridSearchOptions::make_new() {
    vector<double> min = {0.0, 0.0};
    vector<double> max = {10.0, 10.0};
    double expansion = 0.01;
    vector<size_t> ndiv = {5, 5};
    double tol = 1e-4;
    return std::unique_ptr<GridSearchOptions>{new GridSearchOptions{
        min,
        max,
        expansion,
        ndiv,
        tol,
    }};
}

// Allocates new instance
std::unique_ptr<GridSearch> GridSearch::make_new(const std::unique_ptr<GridSearchOptions> &options) {
    // space dimension
    size_t ndim = options->min.size();
    if (ndim < 2 || ndim > 3) {
        throw "min.size() = ndim must be 2 or 3";
    }
    if (options->max.size() != ndim) {
        throw "max.len() must equal ndim = min.len()";
    }

    // expanded borders
    auto vec_min = vector<double>(ndim);
    auto vec_max = vector<double>(ndim);
    for (size_t i = 0; i < ndim; i++) {
        double del = options->max[i] - options->min[i];
        vec_min[i] = options->min[i] - options->expansion * del;
        vec_max[i] = options->max[i] + options->expansion * del;
    }

    // set radius tolerance
    double radius_tol = 0.0;
    for (size_t i = 0; i < ndim; i++) {
        radius_tol += options->tol * options->tol;
    }
    radius_tol = sqrt(radius_tol);

    // compute sizes
    double radius = 0.0;
    auto delta = vector<double>(ndim);
    auto size = vector<double>(ndim);
    for (size_t i = 0; i < ndim; i++) {
        if (options->ndiv[i] < 1) {
            throw "ndiv must be ≥ 1";
        }
        delta[i] = vec_max[i] - vec_min[i];
        if (delta[i] <= 0.0) {
            throw "max must be greater than min";
        }
        size[i] = delta[i] / ((double)options->ndiv[i]);
        if (size[i] <= 2.0 * options->tol) {
            throw "container size = (max-min)/ndiv must be > 2·tol; reduce ndiv or tol";
        }
        radius += size[i] * size[i] / 4.0;
    }
    radius = sqrt(radius);

    // other data
    size_t ncorner = pow(2, ndim);
    auto cf = vector<size_t>(3);  // must be 3
    auto tol = vector<double>(ndim);
    auto halo = vector<vector<double>>(ncorner);
    auto containers = map<Index, map<ID, Item>>();
    for (size_t i = 0; i < ndim; i++) {
        tol[i] = options->tol;
    }
    for (size_t m = 0; m < ncorner; m++) {
        halo[m] = vector<double>(ndim);
    }

    return std::unique_ptr<GridSearch>{
        new GridSearch{
            ndim,
            options->ndiv,
            vec_min,
            vec_max,
            delta,
            size,
            cf,
            tol,
            halo,
            ncorner,
            containers,
        }};
};

// Inserts a new item to the right container in the grid
//
// # Input
//
// * `id` -- identification number for the item
// * `x` -- coordinates (ndim) of the item
void GridSearch::insert(size_t id, vector<double> &x) {
    // check
    if (x.size() != this->ndim) {
        throw "x.size() must equal ndim";
    }

    // add point to container
    int index = this->container_index(x);
    if (index < 0) {
        throw "cannot insert point outside the grid";
    }
    this->update_or_insert(index, id, x);

    // add point to containers touched by halo corners
    this->set_halo(x);
    for (size_t c = 0; c < this->ncorner; c++) {
        int index_corner = this->container_index(this->halo[c]);
        if (index_corner >= 0) {
            if (index_corner != index) {
                this->update_or_insert(index_corner, id, x);  // make sure to use original `x`
            }
        }
    }
}

// Calculates the container index where the point x should be located
//
// # Output
//
// * returns the index of the container or -1 if the point is out-of-range
int GridSearch::container_index(vector<double> &x) {
    auto ratio = vector<size_t>(this->ndim);  // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)
    size_t index = 0;
    for (size_t i = 0; i < this->ndim; i++) {
        if (x[i] < this->min[i] || x[i] > this->max[i]) {
            return -1;  // out-of-range
        }
        ratio[i] = (size_t)((x[i] - this->min[i]) / this->size[i]);
        if (ratio[i] == this->ndiv[i]) {
            // the point is exactly on the max edge, thus select inner container
            ratio[i] -= 1;  // move to the inside
        }
        index += ratio[i] * this->cf[i];
    }
    return index;
}

// Updates container or inserts point in an existing container
void GridSearch::update_or_insert(Index index, ID id, vector<double> &x) {
    auto iter = this->containers.find(index);
    if (iter == this->containers.end()) {
        map<Index, Item> container = {{index, Item{id, x}}};
    } else {
        iter->second.insert({index, Item{id, x}});
    }
}

/// Sets square/cubic halo around point
void GridSearch::set_halo(vector<double> &x) {
    if (this->ndim == 2) {
        this->halo[0][0] = x[0] - this->tol[0];
        this->halo[0][1] = x[1] - this->tol[1];

        this->halo[1][0] = x[0] + this->tol[0];
        this->halo[1][1] = x[1] - this->tol[1];

        this->halo[2][0] = x[0] + this->tol[0];
        this->halo[2][1] = x[1] + this->tol[1];

        this->halo[3][0] = x[0] - this->tol[0];
        this->halo[3][1] = x[1] + this->tol[1];
    } else {
        this->halo[0][0] = x[0] - this->tol[0];
        this->halo[0][1] = x[1] - this->tol[1];
        this->halo[0][2] = x[2] - this->tol[2];

        this->halo[1][0] = x[0] + this->tol[0];
        this->halo[1][1] = x[1] - this->tol[1];
        this->halo[1][2] = x[2] - this->tol[2];

        this->halo[2][0] = x[0] + this->tol[0];
        this->halo[2][1] = x[1] + this->tol[1];
        this->halo[2][2] = x[2] - this->tol[2];

        this->halo[3][0] = x[0] - this->tol[0];
        this->halo[3][1] = x[1] + this->tol[1];
        this->halo[3][2] = x[2] - this->tol[2];

        this->halo[4][0] = x[0] - this->tol[0];
        this->halo[4][1] = x[1] - this->tol[1];
        this->halo[4][2] = x[2] + this->tol[2];

        this->halo[5][0] = x[0] + this->tol[0];
        this->halo[5][1] = x[1] - this->tol[1];
        this->halo[5][2] = x[2] + this->tol[2];

        this->halo[6][0] = x[0] + this->tol[0];
        this->halo[6][1] = x[1] + this->tol[1];
        this->halo[6][2] = x[2] + this->tol[2];

        this->halo[7][0] = x[0] - this->tol[0];
        this->halo[7][1] = x[1] + this->tol[1];
        this->halo[7][2] = x[2] + this->tol[2];
    }
}

void GridSearch::print_details() {
    for (const auto &[index, container] : this->containers) {
        cout << index << ":";
        for (const auto &[id, item] : container) {
            cout << id << "(";
            for (size_t dim = 0; dim < this->ndim; dim++) {
                cout << item.x[dim];
                if (dim < this->ndim - 1) {
                    cout << ",";
                }
            }
            cout << ") ";
        }
        cout << endl;
    }
}