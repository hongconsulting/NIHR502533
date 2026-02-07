// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_SURVFIT_H
#define WH_SURVFIT_H

struct WH_survfit {
  Eigen::VectorXd time;
  Eigen::VectorXd risk;
  Eigen::VectorXd riskESS;
  Eigen::VectorXd event;
  Eigen::VectorXd censor;
  Eigen::VectorXd surv;
  Eigen::VectorXd SE;
  WH_survfit() : time(), risk(), riskESS(), event(), censor(), surv(), SE() {}
};

#endif