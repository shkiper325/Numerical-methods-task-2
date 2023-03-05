#ifndef RUNGE_HPP
#define RUNGE_HPP

#include <vector>
#include <functional>

#include <quadmath.h>

#include "Point.hpp"
#include "number_type.hpp"

struct ApproximatedFunction{
    std::vector<NT> points;
    std::vector<Point> vals;
};

ApproximatedFunction runge5_constant_step(
                                          std::function<Point(NT,Point)> f,
                                          NT a,
                                          NT b,
                                          Point start_val, 
                                          int n
                                         );

ApproximatedFunction runge5_variable_step(
                                          std::function<Point(NT,Point)> f,
                                          NT a,
                                          NT b,
                                          Point start_val,
                                          NT start_h,
                                          NT eps
                                         );

void test_runge5_constant_step();
void test_runge5_variable_step();

#endif