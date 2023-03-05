#include <iostream>

#include <quadmath.h>

#include "runge.hpp"
#include "utils.hpp"
#include "newton.hpp"
#include "number_type.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    // test_runge5_variable_step();
    // test_jacobian_and_inverse3x3();

    if (argc != 3) {
        cerr << "Argument count != 2" << endl;
        exit(1);
    }

    NT alpha;
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
        -0.923399448, //0.1,
        -1.175708690, //0.1
        1.846798897,  //0.1,
        0.01,
        1e-25L,
        1e-20L,
        1e-20L,
        step_count
    );

    return 0;
}

    // main_routine(
    //     alpha,
    //     0.1,
    //     0.1,
    //     0.1,
    //     0.01,
    //     1e-20L,
    //     1e-14L,
    //     1e-14L,
    //     step_count
    // );