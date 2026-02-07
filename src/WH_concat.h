// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_CONCAT_H
#define WH_CONCAT_H

Eigen::VectorXd WH_concat(const Eigen::VectorXd& input1, const Eigen::VectorXd& input2) {
  int n1 = input1.size();
  int n2 = input2.size();
  Eigen::VectorXd output(n1 + n2);
  output.head(n1) = input1;
  output.tail(n2) = input2;
  return output;
}

Eigen::VectorXd WH_concat(const Eigen::VectorXd& input1, const Eigen::VectorXd& input2, const Eigen::VectorXd& input3) {
  int n1 = input1.size();
  int n2 = input2.size();
  int n3 = input3.size();
  Eigen::VectorXd output(n1 + n2 + n3);
  output.head(n1) = input1;
  output.segment(n1, n2) = input2;
  output.tail(n3) = input3;
  return output;
}

Eigen::VectorXd WH_concat(const std::vector<Eigen::VectorXd>& input) {
  int tempint = 0;
  for (int i = 0; i < input.size(); ++i) {
    tempint += input[i].size();
  }
  Eigen::VectorXd output(tempint);
  tempint = 0;
  for (int i = 0; i < input.size(); ++i) {
    output.segment(tempint, input[i].size()) = input[i];
    tempint += input[i].size();
  }
  return output;
}

Eigen::VectorXi WH_concat(const std::vector<Eigen::VectorXi>& input) {
  int tempint = 0;
  for (int i = 0; i < input.size(); ++i) {
    tempint += input[i].size();
  }
  Eigen::VectorXi output(tempint);
  tempint = 0;
  for (int i = 0; i < input.size(); ++i) {
    output.segment(tempint, input[i].size()) = input[i];
    tempint += input[i].size();
  }
  return output;
}

#endif