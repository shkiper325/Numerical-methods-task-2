#include <iostream>
#include <vector>

#include <quadmath.h>

#include "runge.hpp"
#include "utils.hpp"
#include "newton.hpp"
#include "number_type.hpp"

using namespace std;

void print_result(vector<NT> result) {
    printf("p1(0): %18.9lf\n", double(result[0]));
    printf("x1(Pi): %18.9lf\n", double(result[1]));
    printf("B0: %18.9lf\n", double(result[2]));

}

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

    double l1 = 0.1;
    double l2 = 0.1;
    double l3 = 0.1;

    for(int i = 0; i < 51; ++i) {
        double alpha = double(i) / 10.;

        auto result = main_routine(
                                   alpha,
                                   l1,
                                   l2,
                                   l3,
                                   0.01,
                                   1e-25L,
                                   1e-20L,
                                   1e-20L,
                                   10
                                  );

        cout << "###################################" << endl;
        printf("alpha: %9.3lf\n", alpha);
        print_result(result);

        l1 = result[3];
        l2 = result[4];
        l3 = result[5];
    }

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