#ifndef RUNGE_HPP
#define RUNGE_HPP

#include <vector>
#include <functional>

#include "Point.hpp"

struct ApproximatedFunction{
    std::vector<double> points;
    std::vector<Point> vals;
};

ApproximatedFunction runge5_constant_step(
                                          std::function<Point(double,Point)> f,
                                          double a,
                                          double b,
                                          Point start_val, 
                                          int n
                                         );

ApproximatedFunction runge5_variable_step(
                                          std::function<Point(double,Point)> f,
                                          double a,
                                          double b,
                                          Point start_val,
                                          double start_h,
                                          double eps
                                         );

void test_runge5_constant_step();
void test_runge5_variable_step();

#endif