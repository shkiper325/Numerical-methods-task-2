#ifndef POINT_HPP
#define POINT_HPP

#include <algorithm>

#include <cmath>

struct Point {
    double x;
    double p;

    Point(double x_ = 0, double p_ = 0) : x(x_), p(p_) {};

    double linf_norm() {
        return std::max(std::fabs(x), std::fabs(p));
    }

    friend Point operator+(Point one, Point other);
    friend Point operator-(Point one);
    friend Point operator-(Point one, Point other);

    friend Point operator*(double coeff, Point point);
    friend Point operator/(Point point, double coeff);
};

#endif