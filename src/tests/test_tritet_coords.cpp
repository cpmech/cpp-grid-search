#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <vector>

#include "../tritet_coords.h"
#include "check.h"
#include "doctest.h"

using namespace std;

TEST_CASE("TriTetCoords") {
    SUBCASE("Triangle") {
        vector<double> xa = {0.0, 0.0};
        vector<double> xb = {1.0, 0.0};
        vector<double> xc = {0.0, 1.0};
        auto zeta = vector<double>(3);

        vector<double> xp = {0.0, 0.0};
        triangle_coords(zeta, xa, xb, xc, xp);
        vector<double> correct = {1.0, 0.0, 0.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_triangle(zeta));

        xp = {1.0, 0.0};
        triangle_coords(zeta, xa, xb, xc, xp);
        correct = {0.0, 1.0, 0.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_triangle(zeta));

        xp = {0.0, 1.0};
        triangle_coords(zeta, xa, xb, xc, xp);
        correct = {0.0, 0.0, 1.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_triangle(zeta));

        xp = {0.5, 0.5};
        triangle_coords(zeta, xa, xb, xc, xp);
        correct = {0.0, 0.5, 0.5};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_triangle(zeta));

        xp = {1.0, 1.0};
        triangle_coords(zeta, xa, xb, xc, xp);
        correct = {-1.0, 1.0, 1.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_triangle(zeta) == false);

        xp = {1e-15, 1e-15};
        triangle_coords(zeta, xa, xb, xc, xp);
        CHECK(in_triangle(zeta));

        xp = {1e-15, -1e-15};
        triangle_coords(zeta, xa, xb, xc, xp);
        CHECK(in_triangle(zeta) == false);
    }

    SUBCASE("Tetrahedron") {
        vector<double> xa = {0.0, 0.0, 0.0};
        vector<double> xb = {1.0, 0.0, 0.0};
        vector<double> xc = {0.0, 1.0, 0.0};
        vector<double> xd = {0.0, 0.0, 1.0};
        auto zeta = vector<double>(4);

        vector<double> xp = {0.0, 0.0, 0.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        vector<double> correct = {1.0, 0.0, 0.0, 0.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta));

        xp = {1.0, 0.0, 0.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {0.0, 1.0, 0.0, 0.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta));

        xp = {0.0, 1.0, 0.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {0.0, 0.0, 1.0, 0.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta));

        xp = {0.0, 0.0, 1.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {0.0, 0.0, 0.0, 1.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta));

        xp = {0.5, 0.5, 0.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {0.0, 0.5, 0.5, 0.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta));

        xp = {0.5, 0.5, 0.5};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {-0.5, 0.5, 0.5, 0.5};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta) == false);

        xp = {1.0, 1.0, 1.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {-2.0, 1.0, 1.0, 1.0};
        CHECK(equal_vectors(zeta, correct));
        CHECK(in_tetrahedron(zeta) == false);

        xp = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        tetrahedron_coords(zeta, xa, xb, xc, xd, xp);
        correct = {0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        CHECK(equal_vectors_tol(zeta, correct, 1e-15));
        CHECK(in_tetrahedron(zeta));
    }
}