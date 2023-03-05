#include <quadmath.h>

#include "Point.hpp"
#include "number_type.hpp"

Point operator+(Point one, Point other) {
    return Point(one.x + other.x, one.p + other.p);
}

Point operator-(Point one) {
    return Point(-one.x, -one.p);
}

Point operator-(Point one, Point other) {
    return one + (-other);
}

Point operator*(NT coeff, Point point) {
    return Point(coeff * point.x, coeff * point.p);
}

Point operator/(Point point, NT coeff) {
    return Point(point.x / coeff, point.p / coeff);
}