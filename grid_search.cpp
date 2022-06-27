#include "grid_search.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "auxiliary.h"

using namespace std;

// Allocates new instance
std::unique_ptr<GridSearchOptions> GridSearchOptions::make_new() {
    double expansion = 0.01;
    vector<double> xmin = {0.0, 0.0};
    vector<double> xmax = {10.0, 10.0};
    vector<double> tols = {1e-4, 1e-4};
    vector<size_t> ndiv = {5, 5};
    return std::unique_ptr<GridSearchOptions>{new GridSearchOptions{
        expansion,
        xmin,
        xmax,
        tols,
        ndiv,
    }};
}

// Allocates new instance
std::unique_ptr<GridSearch> GridSearch::make_new(const std::unique_ptr<GridSearchOptions> &options) {
    // space dimension
    size_t ndim = options->xmin.size();
    if (ndim < 2 || ndim > 3) {
        throw "min.size() = ndim must be 2 or 3";
    }
    if (options->xmax.size() != ndim) {
        throw "max.len() must equal ndim = min.len()";
    }

    // expanded borders
    auto xmin = vector<double>(ndim);
    auto xmax = vector<double>(ndim);
    for (size_t i = 0; i < ndim; i++) {
        double delta = options->xmax[i] - options->xmin[i];
        xmin[i] = options->xmin[i] - options->expansion * delta;
        xmax[i] = options->xmax[i] + options->expansion * delta;
    }

    // compute sizes
    auto delta = vector<double>(ndim);
    auto side_length = vector<double>(ndim);
    for (size_t i = 0; i < ndim; i++) {
        if (options->ndiv[i] < 1) {
            throw "ndiv must be ≥ 1";
        }
        delta[i] = xmax[i] - xmin[i];
        if (delta[i] <= 0.0) {
            throw "xmax must be greater than xmin";
        }
        side_length[i] = delta[i] / ((double)options->ndiv[i]);
        if (side_length[i] <= 2.0 * options->tols[i]) {
            throw "container side_length = (xmax-xmin)/ndiv must be > 2·tol; reduce ndiv or tol";
        }
    }

    // coefficients
    vector<size_t> coefficient = {1, options->ndiv[0], options->ndiv[0] * options->ndiv[1]};

    // containers
    auto containers = Containers_t();

    return std::unique_ptr<GridSearch>{
        new GridSearch{
            ndim,
            options->expansion,
            xmin,
            xmax,
            options->tols,
            options->ndiv,
            delta,
            side_length,
            coefficient,
            containers,
        }};
};

// Inserts cells (e.g., triangles or tetrahedra) to the correct containers in the grid
void GridSearch::insert_cells(vector<vector<vector<double>>> &mesh) {
    // auxiliary data
    auto cell_xmin = vector<double>(this->ndim);
    auto cell_xmax = vector<double>(this->ndim);
    auto xcorner = vector<double>(this->ndim);

    // loop over all cells
    for (size_t cell_id = 0; cell_id < mesh.size(); cell_id++) {
        // auxiliary reference to cell
        vector<vector<double>> const &cell = mesh[cell_id];

        // check
        size_t nnode = cell.size();
        if (nnode < 3 || nnode > 4) {
            throw "cell.size() == nnode must be 3 (triangle) or 4 (tetrahedron)";
        }

        // loop over the bounding box of the cell
        compute_bounding_box(this->expansion, cell_xmin, cell_xmax, cell);
        if (this->ndim == 2) {
            // lower left
            xcorner[0] = cell_xmin[0];
            xcorner[1] = cell_xmin[1];
            int key = this->_get_container_key(xcorner);
            if (key < 0) {
                throw "lower left corner is outside the grid!";
            }
            this->_update_or_insert_point(key, cell_id, xcorner);

            // lower right
            xcorner[0] = cell_xmax[0];
            xcorner[1] = cell_xmin[1];
            key = this->_get_container_key(xcorner);
            if (key < 0) {
                throw "lower right corner is outside the grid!";
            }
            this->_update_or_insert_point(key, cell_id, xcorner);

            // upper left
            xcorner[0] = cell_xmin[0];
            xcorner[1] = cell_xmax[1];
            key = this->_get_container_key(xcorner);
            if (key < 0) {
                throw "upper left corner is outside the grid!";
            }
            this->_update_or_insert_point(key, cell_id, xcorner);

            // upper right
            xcorner[0] = cell_xmax[0];
            xcorner[1] = cell_xmax[1];
            key = this->_get_container_key(xcorner);
            if (key < 0) {
                throw "upper right corner is outside the grid!";
            }
            this->_update_or_insert_point(key, cell_id, xcorner);
        } else {
            // TODO
            throw "TODO: 3D version";
        }
    }
}

// Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
int GridSearch::find_cell(vector<double> &x, vector<vector<vector<double>>> &mesh) {
    // check
    if (x.size() != this->ndim) {
        throw "x.size() must equal ndim";
    }

    // compute the key of the container where x should be
    int key = this->_get_container_key(x);
    if (key < 0) {
        throw "cannot perform a search with coordinates outside the limits";
    }

    // extract the container data from the map of containers
    auto iter = this->containers.find(key);
    if (iter == this->containers.end()) {
        return -1;  // there is not container set the key corresponding to x
    }

    // check if point is within triangle/tetrahedron
    auto container = iter->second;
    for (const auto &id : container) {
        auto cell = mesh[id];
        if (this->ndim == 2) {
            if (is_point_in_triangle(x, cell)) {
                return id;
            }
        } else {
            if (is_point_in_tetrahedron(x, cell)) {
                return id;
            }
        }
    }

    return -1;  // not found
}

void GridSearch::print_details() {
    cout << "number of non-empty containers = " << this->containers.size() << endl;
    for (const auto &[key, container] : this->containers) {
        cout << "container # " << key << ": items = [";
        bool first = true;
        for (const auto &id : container) {
            if (!first) {
                cout << ", ";
            }
            cout << id;
            first = false;
        }
        cout << "]" << endl;
    }
}

// Calculates the key of the container where the point x should be located
int GridSearch::_get_container_key(vector<double> &x) {
    auto ratio = vector<size_t>(this->ndim);  // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)
    size_t key = 0;
    for (size_t i = 0; i < this->ndim; i++) {
        if (x[i] < this->xmin[i] || x[i] > this->xmax[i]) {
            return -1;  // out-of-range
        }
        ratio[i] = (size_t)((x[i] - this->xmin[i]) / this->side_length[i]);
        if (ratio[i] == this->ndiv[i]) {
            // the point is exactly on the max edge, thus select inner container
            ratio[i] -= 1;  // move to the inside
        }
        key += ratio[i] * this->coefficient[i];
    }
    return key;
}

// Updates a container or inserts a point into an existing container
void GridSearch::_update_or_insert_point(ContainerKey key, ItemID id, vector<double> &x) {
    auto iter = this->containers.find(key);
    if (iter == this->containers.end()) {
        Container_t container = {{id}};
        this->containers.insert({key, container});
    } else {
        Container_t &container = iter->second;
        container.insert({id});
    }
}
