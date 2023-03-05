#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <functional>

#include <quadmath.h>

#include "Point.hpp"
#include "number_type.hpp"

NT C_dist(std::vector<NT> f, std::vector<NT> g);

std::vector<std::vector<NT> > p2vec(std::vector<Point> vec);

void dump_vec(std::vector<NT> vec, std::string path);

NT integrate(std::vector<NT> t, std::vector<NT> x, std::function<NT(NT)> f);

std::vector<NT> jacobian(std::function<std::vector<NT>(NT, NT, NT)> f, NT x, NT y, NT z, NT eps, NT delta);

std::vector<NT> inverse3x3(std::vector<NT> mat);

void print_mat3x3(std::vector<NT> mat);

void test_jacobian_and_inverse3x3();

std::vector<NT> apply_mat_3x3to3 (std::vector<NT> mat, std::vector<NT> vec);

NT maxq(NT a, NT b);

#endif