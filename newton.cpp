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

double Fedorenko_norm(vector<NT> x, vector<NT> F) {
    int N = F.size();
    vector<NT> F_dot(N);

    F_dot[0] = (F[1] - F[0]) / (x[1] - x[0]);
    F_dot[N - 1] = (F[N - 1] - F[N - 2]) / (x[N - 1] - x[N - 2]);
    for (int i = 1; i < N - 1; ++i) {
        F_dot[i] = (F[i + 1] - F[i - 1]) / (x[i + 1] - x[i - 1]);
    }

    double ret = 0;
    for (int i = 0; i < N; ++i) {
        double new_ret = double(fabsq(F[i] / F_dot[i]));
        ret = ret > new_ret ? ret : new_ret;
    }

    return ret;
}

vector<NT> main_routine(
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

    // printf("Start l1: %12.6lf\n", double(l1));
    // printf("Start l2: %12.6lf\n", double(l2));
    // printf("Start l3: %12.6lf\n", double(l3));

    for (int step = 0; step < step_count; ++step) {
        // cout << "##################################################################" << endl;
        // cout << "Step: " << step << endl;

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

        // printf("\n");
        // printf("l1: %18.9lf\n", double(l1));
        // printf("l2: %18.9lf\n", double(l2));
        // printf("l3: %18.9lf\n", double(l3));

        ApproximatedFunction sol = runge5_variable_step(gen_f(alpha, l1, l2), 0, PI, Point(0, l3), start_h, max_runge_error);
        vector<NT> x = p2vec(sol.vals)[0];
        // printf("Fedorenko norm: %9.3e\n", Fedorenko_norm(sol.points, x));
    }

    ApproximatedFunction sol = runge5_variable_step(gen_f(alpha, l1, l2), 0, PI, Point(0, l3), start_h, max_runge_error);

    double p1 = l3 / 2;
    double x1 = sol.vals.back().x;

    vector<NT> x = p2vec(sol.vals)[0];
    int N = x.size();
    vector<NT> x_dot(N);

    x_dot[0] = (x[1] - x[0]) / (sol.points[1] - sol.points[0]);
    x_dot[N - 1] = (x[N - 1] - x[N - 2]) / (sol.points[N - 1] - sol.points[N - 2]);
    for (int i = 1; i < N - 1; ++i) {
        x_dot[i] = (x[i + 1] - x[i - 1]) / (sol.points[i + 1] - sol.points[i - 1]);
    }

    vector<NT> integrant(N);
    for (int i = 0; i < N; ++i) {
        integrant[i] = x_dot[i] * x_dot[i] * double(expq(-alpha * x[i]));
    }

    double B0 = integrate(sol.points, integrant, [](NT t) -> NT {return 1;});

    return {p1, x1, B0, l1, l2, l3};
}