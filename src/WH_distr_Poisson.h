// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_DISTR_POISSON_H
#define WH_DISTR_POISSON_H

#include "WH_RNG.h"

int WH_rPoisson(double lambda) {
  std::poisson_distribution<> dist(lambda);
  return dist(WH_RNG);
}

#endif