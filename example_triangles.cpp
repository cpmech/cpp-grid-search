
#include <cassert>
#include <iostream>

#include "check.h"
#include "grid_search.h"

using namespace std;

int main() {
    try {
        cout << "Running example" << endl;
        cout << "Done" << endl;

    } catch (char const* msg) {
        cout << "ERROR: " << msg << endl;
    } catch (...) {
        cout << "some error occurred" << endl;
    }
    return 0;
}
