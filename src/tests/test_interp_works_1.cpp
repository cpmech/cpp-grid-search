#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../gsearch/grid_search.h"
#include "check.h"

using namespace std;

int main() {
    try {
        // data
        vector<vector<double>> coordinates = {
            {0.0, 0.0, 0.0},  // last column is the temperature
            {0.5, 0.85, 0.986154146165801},
            {1.0, 0.0, 1.0},
            {1.0, 1.7, 1.972308292331602},
            {1.5, 0.85, 1.724093964956667},
            {2.0, 0.0, 2.0},
            {2.0, 1.7, 2.6248809496813372},
            {2.5, 0.85, 2.640549185302179},
            {3.0, 1.7, 3.448187929913334},
        };
        vector<vector<size_t>> triangles = {
            {0, 2, 1},
            {2, 5, 4},
            {1, 2, 4},
            {4, 5, 7},
            {1, 4, 3},
            {4, 7, 6},
            {3, 4, 6},
            {6, 7, 8},
        };

        // allocate grid
        auto grid = GridSearch::make_new(coordinates, triangles);

        // interpolate @ points
        vector<double> x(2);
        for (const auto &x_y_tt : coordinates) {
            x[0] = x_y_tt[0];
            x[1] = x_y_tt[1];
            double tt = grid->find_triangle_and_interpolate(x, coordinates, triangles);
            double error = fabs(tt - x_y_tt[2]);
            assert(error < 1e-17);
        }

        // interpolate @ (1.5,1.0)
        x = {1.5, 1.0};
        double tt = grid->find_triangle_and_interpolate(x, coordinates, triangles);
        double correct = sqrt(x[0] * x[0] + x[1] * x[1]);
        double error = fabs(tt - correct);
        cout << "T = " << tt << " (" << correct << ") error = " << error << endl;
        assert(error < 0.025);

        // handle NaN
        x = {0.5, 1.5};
        tt = grid->find_triangle_and_interpolate(x, coordinates, triangles);
        cout << "T = " << tt << endl;
        assert(isnan(tt));

        cout << "OK" << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
