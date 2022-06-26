#include <cassert>
#include <iostream>

#include "check.h"
#include "grid_search.h"

using namespace std;

int main() {
    try {
        cout << "Running example" << endl;

        auto options = GridSearchOptions::make_new();
        options->min = {0.0, 0.0};
        options->max = {3.0, 3.0};
        options->ndiv = {3, 3};
        auto grid = GridSearch::make_new(options);

        vector<size_t> ndiv_correct = {3, 3};
        assert(equal_vectors(grid->ndiv, ndiv_correct));

        vector<double> delta_correct = {3.06, 3.06};
        assert(equal_vectors_tol(grid->delta, delta_correct, 1e-15));

        vector<double> size_correct = {1.02, 1.02};
        assert(equal_vectors_tol(grid->size, size_correct, 1e-15));

        vector<double> x = {1.0, 2.0};
        grid->insert(123, x);

        grid->print_details();

        cout << "Done" << endl;

    } catch (char const* msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
