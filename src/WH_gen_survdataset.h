// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_GEN_SURVDATASET_H
#define WH_GEN_SURVDATASET_H

#include "WH_gen_recdate_Poisson.h"
#include "WH_gen_survoutcome.h"
#include "WH_randomization.h"
#include "WH_survdata.h"
// #include "WH_survdata_Rcpp.h"

WH_survdata WH_gen_survdataset(Eigen::VectorXd lambda, double lambda_rec,
                               Eigen::VectorXd lambda_LTFU,
                               int n_max,
                               Eigen::VectorXi rand_block, Eigen::VectorXd rand_weight,
                               double t_recmax, double t_max) {
  Eigen::VectorXd recdate = WH_gen_recdate_Poisson(n_max, lambda_rec, t_recmax);
  Eigen::VectorXi arm;
  if (rand_weight(0) != 0.0) { //set randweight to 0 to choose autoblock
    arm = WH_weightrand(rand_weight, n_max);
  } else {
    arm = WH_blockrand(rand_block, n_max);
  }
  Eigen::MatrixXd survoutcome = WH_gen_survoutcome(arm, lambda, lambda_LTFU, t_max);
  Eigen::VectorXd survdate = survoutcome.col(0) + recdate; //column 0 is time, column 1 is status
  Eigen::VectorXd survstatus = survoutcome.col(1);
  WH_survdata output;
  output.recdate = recdate;
  output.arm = arm;
  output.survdate = survdate;
  output.survstatus = survstatus;
  output.survtime = survoutcome.col(0);
  return output;
}

#endif
