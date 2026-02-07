// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_DISTR_BERNOULLI_H
#define WH_DISTR_BERNOULLI_H

#include "WH_RNG.h"

Eigen::VectorXd WH_rBernoulli(int n, double prob) {
  //if (prob < 0.0 || prob > 1.0) throw std::invalid_argument("invalid probability");
  Eigen::VectorXd output(n);
  std::bernoulli_distribution dist(prob);
  for (int i = 0; i < n; ++i) output[i] = dist(WH_RNG);
  return output;
}

Eigen::VectorXd WH_rBernoulli(const Eigen::VectorXd& probs) {
  Eigen::VectorXd output(probs.size());
  for (int i = 0; i < probs.size(); ++i) {
    std::bernoulli_distribution dist(probs[i]);
    output[i] = dist(WH_RNG);
  }
  return output;
}

#endif