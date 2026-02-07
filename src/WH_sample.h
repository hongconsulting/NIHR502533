// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_SAMPLE_H
#define WH_SAMPLE_H

#include "WH_RNG.h"

Eigen::VectorXd WH_sample(const Eigen::VectorXd& x, int n) { //sample without replacement
  if (n > x.size()) {throw std::invalid_argument("[WH_sample] n > x.size()");}
  Eigen::VectorXd output = x;
  std::shuffle(output.data(), output.data() + output.size(), WH_RNG);
  return output.head(n);
}

Eigen::VectorXi WH_sample(const Eigen::VectorXi& x, int n) { //sample without replacement
  if (n > x.size()) {throw std::invalid_argument("[WH_sample] n > x.size()");}
  Eigen::VectorXi output = x;
  std::shuffle(output.data(), output.data() + output.size(), WH_RNG);
  return output.head(n);
}

#endif