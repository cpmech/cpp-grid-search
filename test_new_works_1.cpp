#include <cassert>
#include <iostream>
#include <vector>

#include "grid_search.h"
#include "testing.h"

using namespace std;

int main() {
    try {
        // [num_triangle][nnode=3][ndim=2]
        vector<vector<vector<double>>> tris = {
            {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
            {{1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};

        // lambda function that returns the coordinates of cell's point m
        auto get_x = [&tris](size_t t, size_t m) {
            return &tris[t][m];
        };

        // allocate grid
        auto grid = GridSearch::make_new(tris.size(), get_x);

        // print details
        grid->print_details();

        // check data
        double max_len = 1.0;
        double sl = max_len + 2.0 * GS_TOLERANCE;  // because the bbox is expanded
        assert(grid->side_length == sl);
        vector<size_t> ndiv_correct = {2, 2};
        assert(equal_vectors(grid->ndiv, ndiv_correct));
        vector<double> xmin_correct = {-0.1, -0.1};
        assert(equal_vectors(grid->xmin, xmin_correct));

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
