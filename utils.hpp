#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>

#include "Point.hpp"

long double C_dist(std::vector<long double> f, std::vector<long double> g);

std::vector<std::vector<long double> > p2vec(std::vector<Point> vec);

void dump_vec(std::vector<long double> vec, std::string path);

long double integrate(std::vector<long double> t, std::vector<long double> x, std::function<long double(long double)> f);

std::vector<long double> jacobian(std::function<std::vector<long double>(long double, long double, long double)> f, long double x, long double y, long double z, long double eps, long double delta);

std::vector<long double> inverse3x3(std::vector<long double> mat);

void print_mat3x3(std::vector<long double> mat);

void test_jacobian_and_inverse3x3();

std::vector<long double> apply_mat_3x3to3 (std::vector<long double> mat, std::vector<long double> vec);

#endif