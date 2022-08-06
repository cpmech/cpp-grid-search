#include <cassert>
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

        // find triangle given coords
        vector<double> x = {0.4, 0.2};
        assert(grid->find_triangle(x, coordinates, triangles) == 0);
        x = {1.5, 0.5};
        assert(grid->find_triangle(x, coordinates, triangles) == 1);
        x = {2.0, 1.4};
        assert(grid->find_triangle(x, coordinates, triangles) == 5);
        x = {0.1, 0.5};
        assert(grid->find_triangle(x, coordinates, triangles) == -1);
        x = {2.5, 0.2};
        assert(grid->find_triangle(x, coordinates, triangles) == -1);
        x = {10.0, 1.0};
        auto did_throw = false;
        try {
            grid->find_triangle(x, coordinates, triangles);
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
