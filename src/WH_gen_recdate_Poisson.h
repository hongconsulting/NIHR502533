// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_GEN_RECDATE_POISSON_H
#define WH_GEN_RECDATE_POISSON_H

#include "WH_distr_Poisson.h"

Eigen::VectorXd WH_gen_recdate_Poisson(int n_max, double lambda_rec, double t_recmax) {
  Eigen::VectorXd output = Eigen::VectorXd::Constant(n_max, std::numeric_limits<double>::infinity());
  int currentdate = 1;
  int currentpatient = 1;
  int newrecruits = WH_rPoisson(lambda_rec); // number recruited on day 1
  do {
    if (newrecruits == 0) {
      currentdate++;
      newrecruits = WH_rPoisson(lambda_rec); // generate new recruits for the next date
    } else {
      newrecruits--;
      output(currentpatient - 1) = currentdate;
      currentpatient++;
    }
  } while ((currentdate <= t_recmax) && (currentpatient <= n_max));
  return output;
}

#endif