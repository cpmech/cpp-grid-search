#include <iostream>
#include <vector>

#include "delaunay.h"

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

        // allocate Delaunay structure
        auto del = Delaunay::make_new(coordinates);

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
