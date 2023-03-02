#include "runge.hpp"
#include "utils.hpp"
#include "newton.hpp"

int main() {
    // test_runge5_variable_step();
    // test_jacobian_and_inverse3x3();

    main_routine(
        0,
        0.1,
        0.1,
        0.1,
        0.01,
        1e-5,
        10,
        1e-7,
        1e-7
    );

    return 0;
}