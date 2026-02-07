// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_STRING_INT_H
#define WH_STRING_INT_H

#include "WH_string_hyp_to_min.h"

std::string WH_string_int(int value) {
  return WH_string_hyp_to_min(std::to_string(value));
}

#endif