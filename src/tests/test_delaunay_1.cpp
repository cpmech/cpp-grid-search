#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../delaunay.h"
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
    auto del = Delaunay::make_new(cloud, false);

    // check
    auto npoint = del->coordinates.size();
    auto ntriangle = del->triangles.size();
    CHECK(npoint == 3);
    CHECK(ntriangle == 1);
    CHECK(del->coordinates[0][0] == 0.0);
    CHECK(del->coordinates[0][1] == 0.0);
    CHECK(del->coordinates[1][0] == 1.0);
    CHECK(del->coordinates[1][1] == 0.0);
    CHECK(del->coordinates[2][0] == 0.0);
    CHECK(del->coordinates[2][1] == 1.0);
    CHECK(del->triangles[0][0] == 0);
    CHECK(del->triangles[0][1] == 1);
    CHECK(del->triangles[0][2] == 2);
}
