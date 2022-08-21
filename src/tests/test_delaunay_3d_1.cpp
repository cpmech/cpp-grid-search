#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../delaunay_3d.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("Delaunay") {
    // "cloud" of points
    vector<vector<double>> cloud = {
        {0.0, 0.0, 0.0},  // 0
        {1.0, 0.0, 0.0},  // 1
        {1.0, 1.0, 0.0},  // 2
        {0.0, 1.0, 0.0},  // 3
        {0.0, 0.0, 1.0},  // 4
        {1.0, 0.0, 1.0},  // 5
        {1.0, 1.0, 1.0},  // 6
        {0.0, 1.0, 1.0},  // 7
    };

    // generate Delaunay triangulation
    auto tets = delaunay_3d(cloud, false);

    // check
    auto ntet = tets.size();
    CHECK(ntet == 6);
}
