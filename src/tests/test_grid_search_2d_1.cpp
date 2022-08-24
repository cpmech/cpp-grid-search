#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../grid_search.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("Find") {
    // data
    vector<vector<double>> points = {
        {0.0307942, 0.459123, 1.0},   // 0 last column is the temperature
        {0.0980015, 0.981755, 1.0},   // 1
        {0.133721, 0.348832, 1.0},    // 2
        {0.13928, 0.180603, 1.0},     // 3
        {0.230951, 0.558482, 1.0},    // 4
        {0.478554, 0.00869692, 1.0},  // 5
        {0.540745, 0.331184, 1.0},    // 6
        {0.578587, 0.760349, 1.0},    // 7
        {0.648071, 0.369534, 1.0},    // 8
        {0.903726, 0.975904, 1.0},    // 9
    };
    vector<vector<size_t>> cells = {
        {4, 2, 6},  //  0
        {3, 2, 0},  //  1
        {0, 4, 1},  //  2 << large on y
        {4, 0, 2},  //  3
        {1, 4, 7},  //  4
        {2, 3, 6},  //  5
        {6, 7, 4},  //  6
        {6, 5, 8},  //  7
        {7, 8, 9},  //  8 << very large
        {8, 7, 6},  //  9
        {7, 9, 1},  // 10 << very large
        {6, 3, 5},  // 11
    };

    // allocate grid
    auto grid = GridSearch::make_new(2, points, cells);
    grid->print_details();

    // check initialization data
    double max_len = points[1][1] - points[0][1];
    double sl = max_len + 2.0 * GS_TOLERANCE;  // because the bbox is expanded
    double g = GS_BORDER_TOL;
    CHECK(grid->side_length == sl);
    CHECK(grid->ndiv[0] == 2);
    CHECK(grid->ndiv[1] == 2);
    CHECK(grid->xmin[0] == points[0][0] - g);
    CHECK(grid->xmin[1] == points[5][1] - g);
    CHECK(grid->xmax[0] == points[0][0] - g + sl * 2.0);
    CHECK(grid->xmax[1] == points[5][1] - g + sl * 2.0);
    CHECK(grid->bounding_boxes.size() == cells.size());
    CHECK(grid->containers.size() == 4);
    CHECK(grid->large_cells.size() == 2);
    auto bbox_0 = grid->bounding_boxes[0];
    CHECK(bbox_0[0][0] == points[2][0]);
    CHECK(bbox_0[0][1] == points[6][0]);
    CHECK(bbox_0[1][0] == points[6][1]);
    CHECK(bbox_0[1][1] == points[4][1]);

    // find triangle given coords
    size_t status;
    auto zeta = vector<double>(3);
    vector<double> x = {0.4, 0.2};

    x = {0.4, 0.2};
    auto id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 11);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 1.0, 1e-15));

    x = {0.6, 0.3};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 7);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 1.0, 1e-15));

    x = {0.1, 0.7};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 2);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 1.0, 1e-15));

    x = {0.8, 0.8};
    id = grid->find(status, zeta, x);
    CHECK(status == SUCCESS);
    CHECK(id == 8);
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 1.0, 1e-15));

    id = grid->find(status, zeta, points[5]);
    CHECK(status == SUCCESS);
    if (id != 7) {
        CHECK(id == 11);
    }
    CHECK(equal_scalars_tol(grid->interp(id, zeta), 1.0, 1e-15));

    x = {0.1, 0.1};
    id = grid->find(status, zeta, x);
    CHECK(status == CELL_NOT_FOUND);
    CHECK(id == -1);

    x = {0.6, 0.2};
    id = grid->find(status, zeta, x);
    CHECK(status == CELL_NOT_FOUND);
    CHECK(id == -1);

    x = {0.4, 1.0};
    id = grid->find(status, zeta, x);
    CHECK(status == CELL_NOT_FOUND);
    CHECK(id == -1);

    x = {10.0, 1.0};
    id = grid->find(status, zeta, x);
    CHECK(status == FAIL_POINT_OUTSIDE);
    CHECK(id == -1);
}
