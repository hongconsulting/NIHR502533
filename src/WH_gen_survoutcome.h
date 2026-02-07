// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_GEN_SURVOUTCOME_H
#define WH_GEN_SURVOUTCOME_H

#include "WH_distr_exp.h"

Eigen::MatrixXd WH_gen_survoutcome(const Eigen::VectorXi& arm, const Eigen::VectorXd& lambda, const Eigen::VectorXd& lambda_LTFU, double t_max) {
  int n = arm.size();
  Eigen::MatrixXd output(n, 2);
  for (int i = 0; i < n; i++) {
    double time = WH_rexp(1.0 / lambda[arm[i]]); // parameterised as scale
    double time_LTFU;
    if (lambda_LTFU[arm[i]] > 0.0) {
      time_LTFU = WH_rexp(1.0 / lambda_LTFU[arm[i]]); // parameterised as scale
    } else {
      time_LTFU = t_max + 1; // if 0 LTFU hazard, will never LTFU
    }
    if (time <= time_LTFU) {
      output(i, 0) = time; // column 0 is time
      output(i, 1) = 1; // column 1 is status
    } else {
      output(i, 0) = time_LTFU; // column 0 is time
      output(i, 1) = 0; // column 1 is status
    }
  }
  for (int i = 0; i < n; i++) {
    output(i, 0) = std::ceil(output(i, 0)); // round up
  }
  return output;
}

#endif