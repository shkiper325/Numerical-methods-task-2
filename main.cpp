#include <iostream>

#include "runge.hpp"
#include "utils.hpp"
#include "newton.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    // test_runge5_variable_step();
    // test_jacobian_and_inverse3x3();

    if (argc != 3) {
        cerr << "Argument count != 2" << endl;
        exit(1);
    }

    double alpha;
    int step_count;
    try {
        alpha = atof(argv[1]);
        step_count = atoi(argv[2]);
    }
    catch (...) {
        cerr << "Cant read cmd args" << endl;
        exit(1);
    }

    main_routine(
        alpha,
        0.1,
        0.1,
        0.1,
        0.01,
        1e-5,
        10,
        1e-7,
        1e-7,
        step_count
    );

    return 0;
}