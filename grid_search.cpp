#include "grid_search.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "auxiliary.h"

using namespace std;

// GridSearch tolerance for all directions
const double GS_TOLERANCE = 1e-4;

// GridSearch border tolerance to handle imprecision near the borders
const double GS_BORDER_TOL = 1e-2;

// Allocates new instance
std::unique_ptr<GridSearch> GridSearch::make_new(
    size_t ndim,
    size_t ncell,
    function<size_t(size_t)> get_nnode,
    function<vector<double> const *(size_t, size_t)> get_x) {
    // check input
    if (ndim < 2 || ndim > 3) {
        throw "ndim must be 2 or 3";
    }

    // allocate some variables
    auto xmin = vector<double>(ndim);
    auto xmax = vector<double>(ndim);
    auto x_min_max = BboxMinMax(ndim);
    auto bbox_large = vector<double>(ndim);
    vector<BboxMinMax> bounding_boxes;
    for (size_t i = 0; i < ndim; i++) {
        x_min_max[i] = vector<double>(N_MIN_MAX);
        bbox_large[i] = numeric_limits<double>::min();
    }

    // find limits, bounding boxes, and largest cell
    for (CellId cell_id = 0; cell_id < ncell; cell_id++) {
        size_t nnode = get_nnode(cell_id);
        for (size_t m = 0; m < nnode; m++) {
            auto x = get_x(cell_id, m);
            if (x->size() != ndim) {
                throw "x.size() must be equal to ndim";
            }
            for (size_t i = 0; i < ndim; i++) {
                // limits
                xmin[i] = min(xmin[i], (*x)[i]);
                xmax[i] = max(xmax[i], (*x)[i]);
                // bounding box
                if (m == 0) {
                    for (size_t i = 0; i < ndim; i++) {
                        x_min_max[i][I_MIN] = (*x)[i];
                        x_min_max[i][I_MAX] = (*x)[i];
                    }
                } else {
                    x_min_max[i][I_MIN] = min(x_min_max[i][I_MIN], (*x)[i]);
                    x_min_max[i][I_MAX] = max(x_min_max[i][I_MAX], (*x)[i]);
                }
            }
        }
        // largest cell
        for (size_t i = 0; i < ndim; i++) {
            bbox_large[i] = max(bbox_large[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
        }
        // add to bounding box maps
        bounding_boxes.push_back(x_min_max);
    }

    // make the side_length equal to the largest bounding box dimension
    double side_length = bbox_large[0];
    for (size_t i = 1; i < ndim; i++) {
        side_length = max(side_length, bbox_large[i]);
    }

    // expand side_length by two times the tolerance (left and right)
    side_length += 2.0 * GS_TOLERANCE;

    // expand borders
    for (size_t i = 0; i < ndim; i++) {
        xmin[i] -= GS_BORDER_TOL;
        xmax[i] += GS_BORDER_TOL;
    }

    // number of divisions
    vector<size_t> ndiv(ndim);
    for (size_t i = 0; i < ndim; i++) {
        ndiv[i] = ceil((xmax[i] - xmin[i]) / side_length);
    }

    // update xmax after deciding on the side_length and number of divisions
    for (size_t i = 0; i < ndim; i++) {
        xmax[i] = xmin[i] + side_length * ((double)ndiv[i]);
    }

    // insert cells
    Containers_t containers;
    vector<double> x(ndim);
    for (size_t cell_id = 0; cell_id < ncell; cell_id++) {
        auto x_min_max = bounding_boxes[cell_id];
        for (size_t r = 0; r < 2; r++) {
            for (size_t s = 0; s < 2; s++) {
                for (size_t t = 0; t < (ndim - 1); t++) {
                    x[0] = x_min_max[0][r];
                    x[1] = x_min_max[1][s];
                    if (ndim == 3) {
                        x[2] = x_min_max[2][t];
                    }
                    auto key = calc_container_key(ndim, side_length, ndiv, xmin, x);
                    auto iter = containers.find(key);
                    if (iter == containers.end()) {
                        Container_t container = {{cell_id}};
                        containers.insert({key, container});
                    } else {
                        Container_t &container = iter->second;
                        container.insert({cell_id});
                    }
                }
            }
        }
    }

    // allocate grid
    return std::unique_ptr<GridSearch>{
        new GridSearch{
            ndim,
            ndiv,
            xmin,
            xmax,
            side_length,
            bounding_boxes,
            containers,
        }};
};

// Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
//
// * Returns the CellId or -1 if no cell contains the point
// * `is_in_cell` is a function of cell_id and x that tells whether the point si in the cell or not
int GridSearch::find_cell(vector<double> &x, function<bool(size_t, vector<double> const *)> is_in_cell) {
    // check if the point is out-of-bounds
    for (size_t i = 0; i < this->ndim; i++) {
        if (x[i] < this->xmin[i] || x[i] > this->xmax[i]) {
            throw "given point coordinates are outsize the grid";
        }
    }

    // get the container where `x` falls in
    auto key = calc_container_key(this->ndim, this->side_length, this->ndiv, this->xmin, x);
    auto iter = this->containers.find(key);
    if (iter == this->containers.end()) {
        return -1;  // there is not container set the key corresponding to x
    }

    // find the cell where the point falls in
    auto container = iter->second;
    for (const auto &cell_id : container) {
        auto x_min_max = this->bounding_boxes[cell_id];
        for (size_t i = 0; i < this->ndim; i++) {
            if (x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX]) {
                continue;  // outside the bounding box
            }
        }
        if ((is_in_cell)(cell_id, &x)) {
            return cell_id;
        }
    }

    // not found
    return -1;
}

// Print details about the grid (e.g., for debugging)
void GridSearch::print_details() {
    cout << "number of non-empty containers = " << this->containers.size() << endl;
    for (const auto &[key, container] : this->containers) {
        cout << "container # " << key << ": cells = [";
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
