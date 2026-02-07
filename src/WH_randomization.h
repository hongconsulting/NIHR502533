// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_RANDOMIZATION_H
#define WH_RANDOMIZATION_H

#include "WH_replicate.h"
#include "WH_RNG.h"
#include "WH_sample.h"

Eigen::VectorXi WH_blockrand(const Eigen::VectorXi& block, int n) {
  Eigen::VectorXi answer(n);
  answer.setConstant(-1);
  Eigen::VectorXi block2 = WH_replicate(block, 2);
  int currentpatient = 1;
  int currentblockindex = 0;
  Eigen::VectorXi currentblock;
  do {
    if (currentblockindex == 0) {
      if (WH_rstdunif() < 0.5) {
        currentblock = WH_sample(block, block.size());
        currentblockindex = block.size();
      } else {
        currentblock = WH_sample(block2, block2.size());
        currentblockindex = block2.size();
      }
    }
    answer(currentpatient - 1) = currentblock(currentblockindex - 1);
    currentblockindex--;
    currentpatient++;
  } while (currentpatient <= n);
  return answer;
}

Eigen::VectorXi WH_weightrand(const Eigen::VectorXd& weight, int N) {
  std::discrete_distribution<> dist(weight.data(), weight.data() + weight.size());
  Eigen::VectorXi output(N);
  for (int i = 0; i < N; ++i) {
    output(i) = dist(WH_RNG); //arms start from 0
  }
  return output;
}

#endif