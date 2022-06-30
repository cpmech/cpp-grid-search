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
            {{0.230951, 0.558482}, {0.133721, 0.348832}, {0.540745, 0.331184}},
            {{0.13928, 0.180603}, {0.133721, 0.348832}, {0.0307942, 0.459123}},
            {{0.0307942, 0.459123}, {0.230951, 0.558482}, {0.0980015, 0.981755}},
            {{0.230951, 0.558482}, {0.0307942, 0.459123}, {0.133721, 0.348832}},
            {{0.0980015, 0.981755}, {0.230951, 0.558482}, {0.578587, 0.760349}},
            {{0.133721, 0.348832}, {0.13928, 0.180603}, {0.540745, 0.331184}},
            {{0.540745, 0.331184}, {0.578587, 0.760349}, {0.230951, 0.558482}},
            {{0.540745, 0.331184}, {0.478554, 0.00869692}, {0.648071, 0.369534}},
            {{0.578587, 0.760349}, {0.648071, 0.369534}, {0.903726, 0.975904}},
            {{0.648071, 0.369534}, {0.578587, 0.760349}, {0.540745, 0.331184}},
            {{0.578587, 0.760349}, {0.903726, 0.975904}, {0.0980015, 0.981755}},
            {{0.540745, 0.331184}, {0.13928, 0.180603}, {0.478554, 0.00869692}}};

        // lambda function that returns the coordinates of cell's point m
        auto get_x = [&tris](size_t t, size_t m) {
            return &tris[t][m];
        };

        // allocate grid
        auto grid = GridSearch::make_new(tris.size(), get_x);

        // lambda function that tells whether the point is in the cell or not
        auto is_in_cell = [&tris](size_t t, vector<double> const *x) {
            return in_triangle(tris[t][0], tris[t][1], tris[t][2], (*x));
        };

        // find triangle given coords
        vector<double> x = {0.4, 0.2};
        assert(grid->find_cell(x, is_in_cell) == 11);
        x = {0.6, 0.3};
        assert(grid->find_cell(x, is_in_cell) == 7);
        x = {0.1, 0.7};
        assert(grid->find_cell(x, is_in_cell) == 2);
        x = {0.8, 0.8};
        assert(grid->find_cell(x, is_in_cell) == 8);
        assert(grid->find_cell(tris[7][1], is_in_cell) == 7);
        x = {0.1, 0.1};
        assert(grid->find_cell(x, is_in_cell) == -1);
        x = {0.6, 0.2};
        assert(grid->find_cell(x, is_in_cell) == -1);
        x = {0.4, 1.0};
        assert(grid->find_cell(x, is_in_cell) == -1);
        x = {10.0, 1.0};
        auto did_throw = false;
        try {
            grid->find_cell(x, is_in_cell);
        } catch (const char *msg) {
            did_throw = true;
        }
        assert(did_throw);

        cout << "OK" << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
