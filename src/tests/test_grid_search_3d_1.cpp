#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../grid_search.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("GridSearch 3D - 1") {
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
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 10.0, 1e-15));

    x = {1.0, 0.0, 0.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 20.0, 1e-15));

    x = {0.0, 1.0, 0.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 30.0, 1e-15));

    x = {0.0, 0.0, 1.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 40.0, 1e-15));

    x = {0.5, 0.5, 0.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 25.0, 1e-15));

    x = {0.25, 0.5, 0.25};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 30.0, 1e-15));

    x = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 0);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 30.0, 1e-15));
}
