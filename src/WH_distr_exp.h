// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_DISTR_EXP_H
#define WH_DISTR_EXP_H

#include "WH_RNG.h"

double WH_rexp(double scale) {
  std::exponential_distribution<> dist(1/scale);
  return dist(WH_RNG);
}

double WH_rexp_lambda(double lambda) {
  std::exponential_distribution<> dist(lambda);
  return dist(WH_RNG);
}

Eigen::VectorXd WH_rexp_lambda(int n, double lambda) {
  std::exponential_distribution<> dist(lambda);
  Eigen::VectorXd output(n);
  for (int i = 0; i < n; ++i) output[i] = dist(WH_RNG);
  return output;
}

#endif