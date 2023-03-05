#ifndef POINT_HPP
#define POINT_HPP

#include <algorithm>

#include <cmath>

#include <quadmath.h>

#include "number_type.hpp"

struct Point {
    NT x;
    NT p;

    Point(NT x_ = 0, NT p_ = 0) : x(x_), p(p_) {};

    NT linf_norm() {
        return fabsq(x) > fabsq(p) ? fabsq(x) : fabsq(p);
    }

    friend Point operator+(Point one, Point other);
    friend Point operator-(Point one);
    friend Point operator-(Point one, Point other);

    friend Point operator*(NT coeff, Point point);
    friend Point operator/(Point point, NT coeff);
};

#endif