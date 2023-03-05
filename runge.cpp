#include <functional>
#include <vector>
#include <iostream>
#include <list>

#include <cmath>

#include "runge.hpp"
#include "utils.hpp"

using namespace std;

ApproximatedFunction runge5_constant_step(
                                          std::function<Point(long double, Point)> f,
                                          long double a,
                                          long double b,
                                          Point start_val,
                                          int n)
{
    ApproximatedFunction ret;

    ret.points = vector<long double>(n + 1);
    ret.vals = vector<Point>(n + 1);

    long double h = (b - a) / n;

    ret.points[0] = a;
    ret.vals[0] = start_val;

    for (int i = 0; i < n; ++i) {
        long double xn = a + h * i;
        Point yn = ret.vals[i];

        Point k1 = h * f(xn, yn);
        Point k2 = h * f(xn + h / 2, yn + k1 / 2);
        Point k3 = h * f(xn + h / 2, yn + (k1 + k2) / 4);
        Point k4 = h * f(xn + h, yn - k2 + 2 * k3);
        Point k5 = h * f(xn + 2 * h / 3, yn + (7 * k1 + 10 * k2 + k4));
        Point k6 = h * f(xn + h / 5, yn + (28 * k1 - 125 * k2 + 546 * k3 + 54 * k4 -378 * k5) / 625);

        Point dy = k1 / 24 + 5 * k4 / 48 + 27 * k5 / 56 + 125 * k6 / 336;

        ret.points[i + 1] = xn + h;
        ret.vals[i + 1] = yn + dy;
    }

    return ret;
}

void test_runge5_constant_step() {
    int n = 1000;

    auto f = [](long double t, Point y) -> Point {
        return Point(t + y.p, y.x - t);
    };

    ApproximatedFunction sol = runge5_constant_step(f, 2, 3, Point(1, 2), 1000);

    auto x_target_func = [](long double t) -> long double {return -1 - 3 * exp(2 - t) / 2 + 3 * exp(-2 + t) / 2 + t;};
    auto p_target_func = [](long double t) -> long double {return  1 + 3 * exp(2 - t) / 2 + 3 * exp(-2 + t) / 2 - t;};
    auto x_target = vector<long double>(n + 1);
    auto p_target = vector<long double>(n + 1);
    for (int i = 0; i < n + 1; ++i) {
        x_target[i] = x_target_func(sol.points[i]);
        p_target[i] = p_target_func(sol.points[i]);
    }
   
    auto sol_vals_vectorized = p2vec(sol.vals);
    vector<long double> x_approx = sol_vals_vectorized[0];
    vector<long double> p_approx = sol_vals_vectorized[1];

    cout << "C-dist for x: " << C_dist(x_approx, x_target) << endl;
    cout << "C-dist for p: " << C_dist(p_approx, p_target) << endl;

    dump_vec(x_approx, "tmp/x_approx.txt");
    dump_vec(p_approx, "tmp/p_approx.txt");
}

Point runge_5_dy(
                   std::function<Point(long double,Point)> f,
                   long double t,
                   Point yt,
                   long double h
                  )
{
    Point k1 = h * f(t, yt);
    Point k2 = h * f(t + h / 2, yt + k1 / 2);
    Point k3 = h * f(t + h / 2, yt + (k1 + k2) / 4);
    Point k4 = h * f(t + h, yt - k2 + 2 * k3);
    Point k5 = h * f(t + 2 * h / 3, yt + (7 * k1 + 10 * k2 + k4) / 27);
    Point k6 = h * f(t + h / 5, yt + (28 * k1 - 125 * k2 + 546 * k3 + 54 * k4 - 378 * k5) / 625);

    Point dy = k1 / 24 + (5 * k4) / 48 + (27 * k5) / 56 + (125 * k6) / 336;

    return dy;
}


ApproximatedFunction runge5_variable_step(
                                          std::function<Point(long double,Point)> f,
                                          long double a,
                                          long double b,
                                          Point start_val,
                                          long double start_h,
                                          long double eps
                                         )
{
    auto ret_points = list<long double>();
    auto ret_vals = list<Point>();

    ret_points.push_back(a);
    ret_vals.push_back(start_val);

    long double t = a;
    long double curr_h = start_h;
    Point dy;

    bool working = true;

    int step_counter = 0;

    while (working) {
        if (t + curr_h > b) {
            working = false;
            curr_h = b - t;
        }
        while (true) {
            Point yt = ret_vals.back();

            dy = runge_5_dy(f, t, yt, curr_h);

            auto dy1 = runge_5_dy(f, t, yt, curr_h / 2);
            auto dy2 = runge_5_dy(f, t + curr_h / 2, yt + dy1, curr_h / 2);

            auto error = (dy2 + dy1 - dy).linf_norm();

            if (error > eps) {
                curr_h /= 2;
                working = true;

                continue;
            }
            else if (error < eps / 31) {
                ret_points.push_back(t + curr_h);
                ret_vals.push_back(ret_vals.back() + dy1 + dy2);
                
                t += curr_h;

                curr_h *= 2;

                break;
            }
            else {
                ret_points.push_back(t + curr_h);
                ret_vals.push_back(ret_vals.back() + dy1 + dy2);

                t += curr_h;

                break;
            }
        }

        ++step_counter;
    }

    cout << "Runge-Kutta used " << step_counter << " steps" << endl;

    ApproximatedFunction ret;
    
    ret.points = vector<long double>(ret_points.size());
    ret.vals = vector<Point>(ret_vals.size());
    
    int i = 0;
    for(auto it = ret_points.begin(); it != ret_points.end(); ++it) {
        ret.points[i] = *it;
        ++i;
    }
    i = 0;
    for(auto it = ret_vals.begin(); it != ret_vals.end(); ++it) {
        ret.vals[i] = *it;
        ++i;
    }

    return ret;
}

void test_runge5_variable_step() {
    long double start_h = 0.01;
    long double eps = 1e-5;
    long double K = 3;

    auto f = [](long double t, Point y) -> Point {
        return Point(t + y.p, y.x - t);
    };

    ApproximatedFunction sol = runge5_variable_step(f, 2, 3, Point(1, 2), start_h, eps);

    int n = sol.points.size() - 1;

    auto x_target_func = [](long double t) -> long double {return -1 - 3 * exp(2 - t) / 2 + 3 * exp(-2 + t) / 2 + t;};
    auto p_target_func = [](long double t) -> long double {return  1 + 3 * exp(2 - t) / 2 + 3 * exp(-2 + t) / 2 - t;};
    auto x_target = vector<long double>(n + 1);
    auto p_target = vector<long double>(n + 1);
    for (int i = 0; i < n + 1; ++i) {
        x_target[i] = x_target_func(sol.points[i]);
        p_target[i] = p_target_func(sol.points[i]);
    }
   
    auto sol_vals_vectorized = p2vec(sol.vals);
    vector<long double> x_approx = sol_vals_vectorized[0];
    vector<long double> p_approx = sol_vals_vectorized[1];

    cout << "C-dist for x: " << C_dist(x_approx, x_target) << endl;
    cout << "C-dist for p: " << C_dist(p_approx, p_target) << endl;

    dump_vec(x_approx, "tmp/x_approx.txt");
    dump_vec(p_approx, "tmp/p_approx.txt");
}