// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_DISTR_NORM_H
#define WH_DISTR_NORM_H

#include "WH_constexpr.h"
#include "WH_RNG.h"

double WH_pstdnorm(double x) { // standard normal distribution
  return 0.5 * std::erfc(-x * WH::_SQRT1_2);
}

#endif