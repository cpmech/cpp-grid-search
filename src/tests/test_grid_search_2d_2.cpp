#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../grid_search.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("GridSearch 2D - 2") {
    // data
    vector<vector<double>> points = {
        {0.0, 0.0, 3.0},  // 0 last column is the temperature
        {1.5, 0.0, 2.0},  // 1
        {0.0, 1.0, 1.0},  // 2
    };
    vector<vector<size_t>> cells = {
        {0, 1, 2},  //  0
    };

    // allocate grid
    auto grid = GridSearch::make_new(2, points, cells);
    grid->print_details();

    // output data
    size_t status;
    auto zeta = vector<double>(3);

    // find cell given coords
    vector<double> x = {0.0, 0.0};
    auto id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 3.0, 1e-15));

    x = {1.5, 0.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 2.0, 1e-15));

    x = {0.0, 1.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 1.0, 1e-15));
}
