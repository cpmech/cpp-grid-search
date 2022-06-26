#include <cassert>
#include <iostream>

#include "check.h"
#include "grid_search.h"

using namespace std;

int main() {
    try {
        cout << "Running example" << endl;

        // allocates grid from [0,0] to [3,3] with ndiv = [3,3]
        auto options = GridSearchOptions::make_new();
        options->min = {0.0, 0.0};
        options->max = {3.0, 3.0};
        options->ndiv = {3, 3};
        auto grid = GridSearch::make_new(options);

        // check ndiv
        vector<size_t> ndiv_correct = {3, 3};
        assert(equal_vectors(grid->ndiv, ndiv_correct));

        // check delta. Note the extra 0.06 due to the expansion factor of 1%
        vector<double> delta_correct = {3.06, 3.06};
        assert(equal_vectors_tol(grid->delta, delta_correct, 1e-15));

        // check container/bin size. Note the extra 0.01 due to the expansion factor
        vector<double> size_correct = {1.02, 1.02};
        assert(equal_vectors_tol(grid->size, size_correct, 1e-15));

        // check coefficients [1, ndiv[0], ndiv[0]*ndiv[1]]
        vector<size_t> cf_correct = {1, 3, 9};
        assert(equal_vectors(grid->cf, cf_correct));

        // insert point @ (0.5,0.5) with id = 100
        vector<double> x = {0.5, 0.5};
        grid->insert(100, x);

        // insert point @ (0.8,0.5) with id = 200
        x = {0.8, 0.5};
        grid->insert(200, x);

        // insert point @ (1.0,0.5) with id = 300
        x = {1.0, 0.5};
        grid->insert(300, x);

        // insert point @ (2.9,2.9) with id = 400
        x = {2.9, 2.9};
        grid->insert(400, x);

        // print details
        grid->print_details();

        // find point near (0.5,0.5)
        x = {0.5, 0.5};
        auto id = grid->find(x);
        cout << "found id=" << id << " near (" << x[0] << "," << x[1] << ")" << endl;

        // find point near (1.0,0.5)
        x = {1.0, 0.5};
        id = grid->find(x);
        cout << "found id=" << id << " near (" << x[0] << "," << x[1] << ")" << endl;

        // find point near (2.9,2.9)
        x = {2.9, 2.9};
        id = grid->find(x);
        cout << "found id=" << id << " near (" << x[0] << "," << x[1] << ")" << endl;

        // find point near (3.0,3.0) (=> not found)
        x = {3.0, 3.0};
        id = grid->find(x);
        cout << "found id=" << id << " near (" << x[0] << "," << x[1] << ")" << endl;

        cout << "Done" << endl;

    } catch (char const* msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
