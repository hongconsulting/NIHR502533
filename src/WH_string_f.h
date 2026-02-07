// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_STRING_F_H
#define WH_STRING_F_H

#include "WH_string_hyp_to_min.h"
#include "WH_string_int.h"

std::string WH_string_f(double input, int prec = 2) {
  if (prec < 0) throw std::runtime_error("prec < 0");
  if (prec == 0) return WH_string_int(std::round(input));
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(prec) << input;
  std::string output = oss.str();
  return WH_string_hyp_to_min(output);
}

std::vector<std::string> WH_string_f(const Eigen::VectorXd& input, int prec = 2) {
  std::vector<std::string> output(input.size());
  for (int i = 0; i < input.size(); ++i) {
    output[i] = WH_string_f(input[i], prec);
  }
  return output;
}

std::vector<std::string> WH_string_f(const std::vector<double>& input, int prec = 2) {
  std::vector<std::string> output(input.size());
  for (int i = 0; i < input.size(); ++i) {
    output[i] = WH_string_f(input[i], prec);
  }
  return output;
}

#endif