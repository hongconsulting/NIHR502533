// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_LOGRANK_1SAMPLE_H
#define WH_LOGRANK_1SAMPLE_H

#include "WH_distr_norm.h"

// [[Rcpp::export]]
Eigen::VectorXd WH_expected_survival(const Eigen::VectorXd& survtime, double lambda) {
  return (-lambda * survtime.array()).exp();
}

// [[Rcpp::export]]
double WH_logrank_1sample(Eigen::VectorXd time, Eigen::VectorXd status, double h, bool one_sided = false) {
  Eigen::VectorXd cumhaz = time * h;
  double E = cumhaz.sum();
  double O = status.sum();
  double Z = (O - E) / sqrt(E);
  if (one_sided) {
    return WH_pstdnorm(Z);
  } else {
    return 2 * (1 - WH_pstdnorm(std::abs(Z)));
  }
} // actually the exponential 1-sample score test, used in OneArmTTE

#endif
