// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_DISTR_BINOM_H
#define WH_DISTR_BINOM_H

#include "WH_choose.h"

double WH_dbinom(int k, int n, double prob) {
  return WH_choose(n, k)*std::pow(prob, k)*std::pow(1 - prob, n - k);
}

#endif
