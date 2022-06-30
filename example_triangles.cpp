#include <iostream>
#include <vector>

#include "grid_search.h"

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

        // print details
        grid->print_details();

        // lambda function that tells whether the point is in the cell or not
        auto is_in_cell = [&tris](size_t t, vector<double> const *x) {
            return in_triangle(tris[t][0], tris[t][1], tris[t][2], (*x));
        };

        // find triangle given coords
        vector<double> x = {0.5, 0.5};
        auto id = grid->find_cell(x, is_in_cell);
        cout << "\nx = {" << x[0] << ", " << x[1] << "}" << endl;
        cout << "found triangle with id = " << id << endl;

        // find with another point
        x = {0.4, 0.2};
        id = grid->find_cell(x, is_in_cell);
        cout << "\nx = {" << x[0] << ", " << x[1] << "}" << endl;
        cout << "found triangle with id = " << id << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
