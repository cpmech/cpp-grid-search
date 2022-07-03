#include <cassert>
#include <iostream>
#include <vector>

#include "grid_search.h"
#include "testing.h"

using namespace std;

int main() {
    try {
        // data
        vector<vector<double>> coordinates = {
            {0.0, 0.0},
            {1.0, 0.0},
            {1.0, 1.0},
            {0.0, 1.0},
        };
        vector<vector<size_t>> triangles = {
            {0, 1, 3},
            {1, 2, 3},
        };

        // allocate grid
        auto grid = GridSearch::make_new(coordinates, triangles);

        // check data
        double max_len = 1.0;
        double sl = max_len + 2.0 * GS_TOLERANCE;  // because the bbox is expanded
        assert(grid->side_length == sl);
        vector<size_t> ndiv_correct = {2, 2};
        assert(equal_vectors(grid->ndiv, ndiv_correct));
        vector<double> xmin_correct = {-0.01, -0.01};
        assert(equal_vectors(grid->xmin, xmin_correct));
        vector<double> xmax_correct = {-0.01 + sl * 2.0, -0.01 + sl * 2.0};
        assert(equal_vectors(grid->xmax, xmax_correct));
        assert(grid->bounding_boxes.size() == 2);
        auto bbox_0 = grid->bounding_boxes[0];
        vector<vector<double>> x_min_max_correct = {{0.0, 1.0}, {0.0, 1.0}};
        assert(equal_vectors(bbox_0[0], x_min_max_correct[0]));
        assert(equal_vectors(bbox_0[1], x_min_max_correct[1]));
        auto bbox_1 = grid->bounding_boxes[1];
        assert(equal_vectors(bbox_1[0], x_min_max_correct[0]));
        assert(equal_vectors(bbox_1[1], x_min_max_correct[1]));
        assert(grid->containers.size() == 4);

        cout << "OK" << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
