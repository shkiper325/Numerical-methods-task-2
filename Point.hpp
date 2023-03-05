#ifndef POINT_HPP
#define POINT_HPP

#include <algorithm>

#include <cmath>

struct Point {
    long double x;
    long double p;

    Point(long double x_ = 0, long double p_ = 0) : x(x_), p(p_) {};

    long double linf_norm() {
        return std::max(std::fabs(x), std::fabs(p));
    }

    friend Point operator+(Point one, Point other);
    friend Point operator-(Point one);
    friend Point operator-(Point one, Point other);

    friend Point operator*(long double coeff, Point point);
    friend Point operator/(Point point, long double coeff);
};

#endif