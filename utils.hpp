#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>

#include "Point.hpp"

double C_dist(std::vector<double> f, std::vector<double> g);

std::vector<std::vector<double> > p2vec(std::vector<Point> vec);

void dump_vec(std::vector<double> vec, std::string path);

double integrate(std::vector<double> t, std::vector<double> x, std::function<double(double)> f);

std::vector<double> jacobian(std::function<std::vector<double>(double, double, double)> f, double x, double y, double z, double eps, double delta);

std::vector<double> inverse3x3(std::vector<double> mat);

void print_mat3x3(std::vector<double> mat);

void test_jacobian_and_inverse3x3();

std::vector<double> apply_mat_3x3to3 (std::vector<double> mat, std::vector<double> vec);

#endif