// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_REPLICATE_H
#define WH_REPLICATE_H

Eigen::VectorXd WH_replicate(double input, int times) {
  Eigen::VectorXd output(times);
  output.setConstant(input);  
  return output;
}

Eigen::VectorXd WH_replicate(const Eigen::VectorXd& input, int times) {
  int size = input.size();
  Eigen::VectorXd output(size * times);
  for (int i = 0; i < times; ++i) {
    output.segment(i * size, size) = input;
  }
  return output;
}

Eigen::VectorXi WH_replicate(const Eigen::VectorXi& input, int times) {
  int size = input.size();
  Eigen::VectorXi output(size * times);
  for (int i = 0; i < times; ++i) {
    output.segment(i * size, size) = input;
  }
  return output;
}

#endif