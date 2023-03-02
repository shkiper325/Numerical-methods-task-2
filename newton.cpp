#include <functional>
#include <vector>
#include <iostream>

#include <cmath>

#include "Point.hpp"
#include "runge.hpp"
#include "utils.hpp"

using namespace std;

const double PI = 3.1415;

function<Point(double,Point)> gen_f(double alpha, double l1, double l2) {
    return [alpha = alpha, l1 = l1, l2 = l2](double t, Point y) -> Point {
        double x = y.x;
        double p = y.p;

        Point ret;
        ret.x = p / (2 * exp(-alpha * x));
        ret.p = -alpha * p * p / (4 * exp(-alpha * x)) + l1 * sin(t) + l2 * cos(t) / (1 + alpha * t * t);

        return ret;
    };
}

void main_routine(
                  double alpha,
                  double l1_start,
                  double l2_start,
                  double l3_start,
                  double start_h,
                  double max_runge_error,
                  double K,
                  double eps,
                  double delta
                 )
{
    double l1 = l1_start;
    double l2 = l2_start;
    double l3 = l3_start;

    printf("Start l1: %12.6lf\n", l1);
    printf("Start l2: %12.6lf\n", l2);
    printf("Start l3: %12.6lf\n", l3);


    for (int step = 0; step < 10; ++step) {
        cout << "##################################################################" << endl;
        cout << "Step: " << step << endl;

        auto F = [alpha = alpha, start_h = start_h, max_runge_error = max_runge_error, K = K](double l1, double l2, double l3) -> vector<double> {
            ApproximatedFunction sol = runge5_variable_step(gen_f(alpha, l1, l2), 0, PI, Point(0, l3), start_h, max_runge_error, K);
            
            vector<vector<double> > sol_vectorized = p2vec(sol.vals);

            double F1 = integrate(sol.points, sol_vectorized[0], [](double t) -> double {return sin(t);}) - 1;
            double F2 = integrate(sol.points, sol_vectorized[0], [alpha = alpha](double t) -> double {return cos(t) / (1 + alpha * t * t);});
            double F3 = sol_vectorized[1].back();

            return {F1, F2, F3};
        };

        auto mat = jacobian(F, l1, l2, l3, 1e-7, 1e-7);
        mat = inverse3x3(mat);
        vector<double> minus_dx = apply_mat_3x3to3(mat, F(l1, l2, l3));

        l1 -= minus_dx[0];
        l2 -= minus_dx[1];
        l3 -= minus_dx[2];

        printf("\n");
        printf("l1: %12.6lf\n", l1);
        printf("l2: %12.6lf\n", l2);
        printf("l3: %12.6lf\n", l3);
    }

    printf("\nDone!\n");
}