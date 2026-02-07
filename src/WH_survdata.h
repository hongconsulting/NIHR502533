// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_SURVDATA_H
#define WH_SURVDATA_H

struct WH_survdata {
  Eigen::VectorXd recdate;
  Eigen::VectorXi arm;
  Eigen::VectorXi cluster;
  Eigen::VectorXd survdate;
  Eigen::VectorXd survtime;
  Eigen::VectorXd survstatus;
  Eigen::MatrixXd covar;
  WH_survdata()
    : recdate(), arm(), cluster(), survdate(), survtime(), survstatus(), covar() {}
};

#endif