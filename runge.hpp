#ifndef RUNGE_HPP
#define RUNGE_HPP

#include <vector>
#include <functional>

#include "Point.hpp"

struct ApproximatedFunction{
    std::vector<long double> points;
    std::vector<Point> vals;
};

ApproximatedFunction runge5_constant_step(
                                          std::function<Point(long double,Point)> f,
                                          long double a,
                                          long double b,
                                          Point start_val, 
                                          int n
                                         );

ApproximatedFunction runge5_variable_step(
                                          std::function<Point(long double,Point)> f,
                                          long double a,
                                          long double b,
                                          Point start_val,
                                          long double start_h,
                                          long double eps
                                         );

void test_runge5_constant_step();
void test_runge5_variable_step();

#endif