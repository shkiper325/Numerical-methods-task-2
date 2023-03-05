#include "Point.hpp"

Point operator+(Point one, Point other) {
    return Point(one.x + other.x, one.p + other.p);
}

Point operator-(Point one) {
    return Point(-one.x, -one.p);
}

Point operator-(Point one, Point other) {
    return one + (-other);
}

Point operator*(long double coeff, Point point) {
    return Point(coeff * point.x, coeff * point.p);
}

Point operator/(Point point, long double coeff) {
    return Point(point.x / coeff, point.p / coeff);
}