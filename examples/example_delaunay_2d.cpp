#include <iostream>
#include <vector>

#include "../src/delaunay_2d.h"

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
        auto triangles = delaunay_2d(cloud, false);

        // print the triangles list
        auto ntriangle = triangles.size();
        cout << "ntriangle = " << ntriangle << endl;
        for (size_t t = 0; t < ntriangle; t++) {
            cout << "v0, v1, v2 = " << triangles[t][0] << ", " << triangles[t][1] << ", " << triangles[t][2] << endl;
        }

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
