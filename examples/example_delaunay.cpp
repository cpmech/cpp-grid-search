#include <iostream>
#include <vector>

#include "../src/delaunay.h"

using namespace std;

int main() {
    try {
        // "cloud" of points
        vector<vector<double>> cloud = {
            {0.478554, 0.00869692},
            {0.13928, 0.180603},
            {0.578587, 0.760349},
            {0.903726, 0.975904},
            {0.0980015, 0.981755},
            {0.133721, 0.348832},
            {0.648071, 0.369534},
            {0.230951, 0.558482},
            {0.0307942, 0.459123},
            {0.540745, 0.331184},
        };

        // generate Delaunay triangulation
        auto del = Delaunay::make_new(cloud, false);

        // print the point coordinates
        auto npoint = del->coordinates.size();
        cout << "npoint = " << npoint << endl;
        for (size_t p = 0; p < npoint; p++) {
            cout << "x, y, T = " << del->coordinates[p][0] << ", " << del->coordinates[p][1] << ", " << del->coordinates[p][2] << endl;
        }

        // print the triangle connectivities
        auto ntriangle = del->triangles.size();
        cout << "ntriangle = " << ntriangle << endl;
        for (size_t t = 0; t < ntriangle; t++) {
            cout << "v0, v1, v2 = " << del->triangles[t][0] << ", " << del->triangles[t][1] << ", " << del->triangles[t][2] << endl;
        }

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
