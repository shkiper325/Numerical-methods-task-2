#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>

#include "Point.hpp"

double C_dist(std::vector<double> f, std::vector<double> g);

std::vector<std::vector<double> > p2vec(std::vector<Point> vec);

void dump_vec(std::vector<double> vec, std::string path);

#endif