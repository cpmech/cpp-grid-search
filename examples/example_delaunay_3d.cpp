#include <iostream>
#include <vector>

#include "../src/delaunay_3d.h"

using namespace std;

int main() {
    try {
        // "cloud" of points
        vector<vector<double>> cloud = {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 1.0},
            {0.0, 1.0, 1.0},
        };

        // generate Delaunay tetrahedralization
        auto tets = delaunay_3d(cloud, false);

        // print the tetrahedra list
        auto ntet = tets.size();
        cout << "ntet = " << ntet << endl;
        for (size_t t = 0; t < ntet; t++) {
            cout << "v0, v1, v2, v3 = " << tets[t][0] << ", " << tets[t][1] << ", " << tets[t][2] << ", " << tets[t][3] << endl;
        }

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
