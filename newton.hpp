#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <vector>

#include <quadmath.h>

#include "number_type.hpp"

std::vector<NT> main_routine(
                  NT alpha,
                  NT l1_start,
                  NT l2_start,
                  NT l3_start,
                  NT start_h,
                  NT max_runge_error,
                  NT eps,
                  NT delta,
                  int step_count
                 );

#endif