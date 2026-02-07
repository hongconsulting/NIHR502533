// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_ANALYSE_FIXEDT_SUMMARY_H
#define WH_ANALYSE_FIXEDT_SUMMARY_H

#include "WH_survdata.h"
// #include "WH_survdata_Rcpp.h"
#include "WH_survfit.h"
// #include "WH_survfit_Rcpp.h"

WH_survfit WH_analyse_fixedt_summary(int n_arms, double fixed_time, WH_survdata look) {
  int n = look.arm.size();
  WH_survfit output;
  Eigen::VectorXd FTevent = Eigen::VectorXd::Zero(n_arms);
  Eigen::VectorXd FTrisk = Eigen::VectorXd::Zero(n_arms); // non-events
  Eigen::VectorXd FTriskESS = Eigen::VectorXd::Zero(n_arms);
  if (n > 0) { // check if no one recruited
    for (int i = 0; i < n; i++) {
      if (look.survtime[i] > fixed_time) {
        FTrisk[look.arm[i]]++; // alive after fixed time
        FTriskESS[look.arm[i]]++; // alive after fixed time
      } else if (look.survtime[i] == fixed_time) {
        if (look.survstatus[i]==1) {
          FTevent[look.arm[i]]++; // event at fixed time
        } else {
          FTrisk[look.arm[i]]++; // alive at fixed time
          FTriskESS[look.arm[i]]++; // alive at fixed time
        }
      } else if (look.survtime[i] < fixed_time) {
        if (look.survstatus[i]==1) {
          FTevent[look.arm[i]]++; // event before fixed time
        } else {
          FTriskESS[look.arm[i]] = FTriskESS[look.arm[i]] + look.survtime[i]/fixed_time; // partial observation
        }
      }
    }
  }
  output.event = FTevent;
  output.risk = FTrisk;
  output.riskESS = FTriskESS.array().floor();
  return output; // arm-wise summary statistics for Bayesian analysis, GPC, etc
}

#endif