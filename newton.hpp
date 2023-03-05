#ifndef NEWTON_HPP
#define NEWTON_HPP

void main_routine(
                  long double alpha,
                  long double l1_start,
                  long double l2_start,
                  long double l3_start,
                  long double start_h,
                  long double max_runge_error,
                  long double eps,
                  long double delta,
                  int step_count
                 );

#endif