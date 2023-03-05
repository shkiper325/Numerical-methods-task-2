#include <functional>
#include <vector>
#include <iostream>

#include <cmath>

#include "Point.hpp"
#include "runge.hpp"
#include "utils.hpp"

using namespace std;

const long double PI = 3.1415926535897932384626433832795;

function<Point(long double,Point)> gen_f(long double alpha, long double l1, long double l2) {
    return [alpha = alpha, l1 = l1, l2 = l2](long double t, Point y) -> Point {
        long double x = y.x;
        long double p = y.p;

        Point ret;
        ret.x = p / (2 * exp(-alpha * x));
        ret.p = -alpha * p * p / (4 * exp(-alpha * x)) + l1 * sin(t) + l2 * cos(t) / (1 + alpha * t * t);

        return ret;
    };
}

void main_routine(
                  long double alpha,
                  long double l1_start,
                  long double l2_start,
                  long double l3_start,
                  long double start_h,
                  long double max_runge_error,
                  long double eps,
                  long double delta,
                  int step_count
                 )
{
    long double l1 = l1_start;
    long double l2 = l2_start;
    long double l3 = l3_start;

    printf("Start l1: %12.6Lf\n", l1);
    printf("Start l2: %12.6Lf\n", l2);
    printf("Start l3: %12.6Lf\n", l3);


    for (int step = 0; step < step_count; ++step) {
        cout << "##################################################################" << endl;
        cout << "Step: " << step << endl;

        auto F = [alpha = alpha, start_h = start_h, max_runge_error = max_runge_error](long double l1, long double l2, long double l3) -> vector<long double> {
            ApproximatedFunction sol = runge5_variable_step(gen_f(alpha, l1, l2), 0, PI, Point(0, l3), start_h, max_runge_error);
            
            vector<vector<long double> > sol_vectorized = p2vec(sol.vals);

            long double F1 = integrate(sol.points, sol_vectorized[0], [](long double t) -> long double {return sin(t);}) - 1.;
            long double F2 = integrate(sol.points, sol_vectorized[0], [alpha = alpha](long double t) -> long double {return cos(t) / (1. + alpha * t * t);});
            long double F3 = sol_vectorized[1].back();

            return {F1, F2, F3};
        };

        auto mat = jacobian(F, l1, l2, l3, eps, delta);
        mat = inverse3x3(mat);
        vector<long double> minus_dx = apply_mat_3x3to3(mat, F(l1, l2, l3));

        l1 -= minus_dx[0];
        l2 -= minus_dx[1];
        l3 -= minus_dx[2];

        printf("\n");
        printf("l1: %18.9Lf\n", l1);
        printf("l2: %18.9Lf\n", l2);
        printf("l3: %18.9Lf\n", l3);
    }

    printf("\nDone!\n");
}