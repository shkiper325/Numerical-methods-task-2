#include <functional>
#include <vector>
#include <iostream>

#include <cmath>

#include <quadmath.h>

#include "Point.hpp"
#include "runge.hpp"
#include "utils.hpp"
#include "number_type.hpp"

using namespace std;

const NT PI = 3.1415926535897932384626433832795L;

function<Point(NT,Point)> gen_f(NT alpha, NT l1, NT l2) {
    return [alpha = alpha, l1 = l1, l2 = l2](NT t, Point y) -> Point {
        NT x = y.x;
        NT p = y.p;

        Point ret;
        ret.x = p / (2 * expq(-alpha * x));
        ret.p = -alpha * p * p / (4 * expq(-alpha * x)) + l1 * sinq(t) + l2 * cosq(t) / (1 + alpha * t * t);

        return ret;
    };
}

void main_routine(
                  NT alpha,
                  NT l1_start,
                  NT l2_start,
                  NT l3_start,
                  NT start_h,
                  NT max_runge_error,
                  NT eps,
                  NT delta,
                  int step_count
                 )
{
    NT l1 = l1_start;
    NT l2 = l2_start;
    NT l3 = l3_start;

    printf("Start l1: %12.6lf\n", double(l1));
    printf("Start l2: %12.6lf\n", double(l2));
    printf("Start l3: %12.6lf\n", double(l3));


    for (int step = 0; step < step_count; ++step) {
        cout << "##################################################################" << endl;
        cout << "Step: " << step << endl;

        auto F = [alpha = alpha, start_h = start_h, max_runge_error = max_runge_error](NT l1, NT l2, NT l3) -> vector<NT> {
            ApproximatedFunction sol = runge5_variable_step(gen_f(alpha, l1, l2), 0, PI, Point(0, l3), start_h, max_runge_error);
            
            vector<vector<NT> > sol_vectorized = p2vec(sol.vals);

            NT F1 = integrate(sol.points, sol_vectorized[0], [](NT t) -> NT {return sinq(t);}) - 1.;
            NT F2 = integrate(sol.points, sol_vectorized[0], [alpha = alpha](NT t) -> NT {return cosq(t) / (1. + alpha * t * t);});
            NT F3 = sol_vectorized[1].back();

            return {F1, F2, F3};
        };

        auto mat = jacobian(F, l1, l2, l3, eps, delta);
        mat = inverse3x3(mat);
        vector<NT> minus_dx = apply_mat_3x3to3(mat, F(l1, l2, l3));

        l1 -= minus_dx[0];
        l2 -= minus_dx[1];
        l3 -= minus_dx[2];

        printf("\n");
        printf("l1: %18.9lf\n", double(l1));
        printf("l2: %18.9lf\n", double(l2));
        printf("l3: %18.9lf\n", double(l3));
    }

    printf("\nDone!\n");
}