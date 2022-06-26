#include "grid_search.h"

#include <cmath>
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
    auto containers = map<size_t, map<size_t, Item>>();
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