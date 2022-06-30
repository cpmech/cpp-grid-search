#include <functional>
#include <iostream>
#include <vector>

using namespace std;

const size_t NNODE = 3;

void print_x(size_t num_triangle, function<vector<double> const &(size_t, size_t)> get_x) {
    for (size_t t = 0; t < num_triangle; t++) {
        for (size_t m = 0; m < NNODE; m++) {
            auto x = get_x(t, m);
            cout << "triangle # " << t << ": x" << m << " = ";
            cout << x[0] << "," << x[1] << endl;
        }
    }
}

int main() {
    try {
        // [num_triangle][nnode=3][ndim=2]
        vector<vector<vector<double>>> triangles = {
            {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
            {{1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};

        // lambda function that returns the coordinates of cell's point i
        auto get_x = [&triangles](size_t t, size_t i) {
            return triangles[t][i];
        };

        // print data
        print_x(triangles.size(), get_x);

    } catch (char const *msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
