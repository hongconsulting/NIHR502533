// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_RNG_H
#define WH_RNG_H

#include <random> //RNG

//std::random_device rd;
//thread_local std::mt19937 WH_RNG(rd());
//std::mt19937 WH_RNG(rd());
thread_local std::mt19937 WH_RNG;

inline double WH_rstdunif() {
  std::uniform_real_distribution<> dis(0.0, 1.0);
  return dis(WH_RNG);
}

void WH_set_seed(int seed) {
  //Rcpp::Rcout << "Thread ID: " << std::this_thread::get_id() << std::endl;
  WH_RNG.seed(seed);
}

#endif