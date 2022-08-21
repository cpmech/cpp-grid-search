#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../delaunay_2d.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("Delaunay") {
    // "cloud" of points
    vector<vector<double>> cloud = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
    };

    // generate Delaunay triangulation
    auto triangles = delaunay_2d(cloud, false);

    // check
    auto ntriangle = triangles.size();
    CHECK(ntriangle == 1);
    CHECK(triangles[0][0] == 0);
    CHECK(triangles[0][1] == 1);
    CHECK(triangles[0][2] == 2);
}
