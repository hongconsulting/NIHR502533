// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_ASSERT_H
#define WH_ASSERT_H

bool WH_is_nondecreasing(const Eigen::VectorXd& input) {
  for (int i = 1; i < input.size(); ++i) {
    if (!(input[i] >= input[i - 1])) return false;
  }
  return true;
}

#endif