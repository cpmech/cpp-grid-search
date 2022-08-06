#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>
#include <iostream>
#include <vector>

#include "../gsearch/grid_search.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("Interpolation") {
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

    SUBCASE("find_triangle_and_interpolate works") {
        // interpolate @ points
        vector<double> x(2);
        for (const auto &x_y_tt : coordinates) {
            x[0] = x_y_tt[0];
            x[1] = x_y_tt[1];
            double tt = grid->find_triangle_and_interpolate(x, coordinates, triangles);
            double error = fabs(tt - x_y_tt[2]);
            CHECK(error < 1e-17);
        }

        // interpolate @ (1.5,1.0)
        x = {1.5, 1.0};
        double tt = grid->find_triangle_and_interpolate(x, coordinates, triangles);
        double correct = sqrt(x[0] * x[0] + x[1] * x[1]);
        double error = fabs(tt - correct);
        // cout << "T = " << tt << " (" << correct << ") error = " << error << endl;
        CHECK(error < 0.025);
    }

    SUBCASE("find_triangle_and_interpolate handles NaNs") {
        // handle NaN
        vector<double> x = {0.5, 1.5};
        double tt = grid->find_triangle_and_interpolate(x, coordinates, triangles);
        // cout << "T = " << tt << endl;
        CHECK(isnan(tt));
    }
}
