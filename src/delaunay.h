#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

extern "C" {
#include "triangle/interface_triangle.h"
}

using namespace std;

struct Delaunay {
    ExtTriangle* ext_triangle;  // data allocated by the c-code

    // destructor to make sure we clean Triangle data
    ~Delaunay() {
        std::cout << "Task::Destructor" << std::endl;
    }

    static std::unique_ptr<Delaunay> make_new(vector<vector<double>> const& coordinates) {
        // constants
        size_t npoint = coordinates.size();
        if (npoint < 3) {
            throw "number of points must be >= 3";
        }
        size_t ncol = coordinates[0].size();
        if (ncol < 2) {
            throw "coordinates.ncol must be >= 2";
        }

        // allocate triangle data
        auto ext_triangle = new_triangle(npoint, 0, 0, 0);

        // allocate Delaunay structure
        return std::unique_ptr<Delaunay>{
            new Delaunay{
                ext_triangle,
            }};
    }
};
