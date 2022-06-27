#include <cassert>
#include <iostream>
#include <vector>

#include "auxiliary.h"
#include "grid_search.h"

using namespace std;

int main() {
    try {
        cout << "Running example" << endl;

        // allocates grid from [-1,-1] to [2,2] with ndiv = [10,10]
        auto options = GridSearchOptions::make_new();
        options->xmin = {-1.0, -1.0};
        options->xmax = {2.0, 2.0};
        options->ndiv = {10, 10};
        auto grid = GridSearch::make_new(options);

        // [num_triangle][nnode=3][ndim=2]
        vector<vector<vector<double>>> triangles = {
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

        // add triangles (aka cells) to grid
        grid->insert_cells(triangles);

        // print details
        grid->print_details();

        // find triangle given coords
        vector<double> x = {0.5, 0.5};
        auto id = grid->find_cell(x, triangles);
        cout << "found triangle with id = " << id << endl;

        cout << "Done" << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
