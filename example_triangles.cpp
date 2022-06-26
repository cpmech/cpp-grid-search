#include <cassert>
#include <iostream>
#include <vector>

#include "auxiliary.h"
#include "grid_search.h"

using namespace std;

const size_t NNODE = 3;
const size_t NDIM = 2;

void triangle_bounding_box(vector<double> &xmin, vector<double> &xmax, vector<vector<double>> &triangle) {
    for (size_t i = 0; i < NDIM; i++) {
        xmin[i] = triangle[0][i];
        xmax[i] = triangle[0][i];
        for (size_t m = 1; m < NNODE; m++) {
            xmin[i] = min(xmin[i], triangle[m][i]);
            xmax[i] = max(xmax[i], triangle[m][i]);
        }
    }
}

int main() {
    try {
        cout << "Running example" << endl;

        // allocates grid from [-1,-1] to [2,2] with ndiv = [10,10]
        auto options = GridSearchOptions::make_new();
        options->min = {-1.0, -1.0};
        options->max = {2.0, 2.0};
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

        size_t num_triangle = triangles.size();

        auto x = vector<double>(NDIM);
        auto xmin = vector<double>(NDIM);
        auto xmax = vector<double>(NDIM);

        for (size_t t = 0; t < num_triangle; t++) {
            triangle_bounding_box(xmin, xmax, triangles[t]);
            for (size_t m = 1; m < NNODE; m++) {
                for (size_t i = 0; i < NDIM; i++) {
                    x[i] = triangles[t][m][i];
                }
                grid->insert(t, x);
            }
        }

        // print details
        grid->print_details();

        // find triangle given coords
        x = {0.5, 0.5};
        int ci = grid->container_index(x);
        if (ci < 0) {
            throw "point is outside grid";
        }
        auto iter = grid->containers.find(ci);
        if (iter == grid->containers.end()) {
            throw "cannot find containers with triangles under x";
        }
        auto container = iter->second;
        for (const auto &[t, _] : container) {
            if (is_point_in_triangle(x, triangles[t])) {
                cout << "point is in triangle # " << t << endl;
            }
        }

        cout << "Done" << endl;

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
