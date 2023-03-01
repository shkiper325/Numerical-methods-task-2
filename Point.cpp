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

Point operator*(double coeff, Point point) {
    return Point(coeff * point.x, coeff * point.p);
}

Point operator/(Point point, double coeff) {
    return Point(point.x / coeff, point.p / coeff);
}