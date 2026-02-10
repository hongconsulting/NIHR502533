// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_CHOOSE_H
#define WH_CHOOSE_H

double WH_choose(int n, int k) {
  if (k < 2) {
    if (k <	0) return 0.0;
    if (k == 0) return 1.0;
    return n; // k == 1
  }
  if (n < 0) return WH_choose(- n + k - 1, k);
  if (n < k) return 0.0;
  //
  if (k > n - k) k = n - k; // symmetry
  double result = 1;
  for (int i = 1; i <= k; ++i) {
    result *= (n - i + 1)/static_cast<double>(i); // affects rounding
  }
  return result;
}

#endif
