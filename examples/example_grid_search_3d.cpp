#include <cassert>
#include <iostream>
#include <vector>

#include "../src/grid_search.h"

using namespace std;

int main() {
    try {
        // data
        vector<vector<double>> points = {
            {0.0, 0.0, 0.0, 10.0},  // 0 last column is the temperature
            {1.0, 0.0, 0.0, 20.0},  // 1
            {0.0, 1.0, 0.0, 30.0},  // 2
            {0.0, 0.0, 1.0, 40.0},  // 3
        };
        vector<vector<size_t>> cells = {
            {0, 1, 2, 3},  //  0
        };

        // allocate grid
        size_t ndim = 3;
        auto grid = GridSearch::make_new(ndim, points, cells);
        grid->print_details();

        // output data
        size_t status;
        auto zeta = vector<double>(ndim + 1);

        // find cell given coords
        vector<double> x = {0.0, 0.0, 0.0};
        auto id = grid->find(status, zeta, x);
        assert(status == SUCCESS);
        assert(id == 0);

        // interpolate the temperature @ {0,0,0}
        auto temp = grid->interp(id, zeta);
        cout << "\nx = {" << x[0] << ", " << x[1] << ", " << x[2] << "}" << endl;
        cout << "temperature = " << temp << " (should print 10)" << endl;
        cout << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
