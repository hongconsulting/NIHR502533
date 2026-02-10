// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_BINOM_TEST_H
#define WH_BINOM_TEST_H

#include "WH_distr_binom.h"

double WH_binom_test(int x, int n, double prop_null) {
  double prob = WH_dbinom(x, n, prop_null);
  double output = 0.0;
  for (int k = 0; k <= n; ++k) {
    if (WH_dbinom(k, n, prop_null) <= prob) output += WH_dbinom(k, n, prop_null);
  }
  return std::min(1.0, output);
}

#endif
