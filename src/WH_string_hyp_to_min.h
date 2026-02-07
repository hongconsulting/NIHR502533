// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_STRING_HYP_TO_MIN_H
#define WH_STRING_HYP_TO_MIN_H

std::string WH_string_hyp_to_min(const std::string& input) {
  if (input == "nan") return "NaN";
  std::string output = input;
  const std::string hyphen = "-";
  const std::string minus = "−"; // UTF-8 encoded
  size_t pos = output.find(hyphen);
  while (pos != std::string::npos) {
    output.replace(pos, hyphen.length(), minus);
    pos = output.find(hyphen, pos + minus.length());
  }
  return output;
}

#endif