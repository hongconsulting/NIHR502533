// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_LOOK_SURVDATA_H
#define WH_LOOK_SURVDATA_H

#include "WH_assert.h"
#include "WH_survdata.h"
// #include "WH_survdata_Rcpp.h"

// data look at interim dates
WH_survdata WH_look_survdata(double date, WH_survdata input) { // VectorXd version
  if (!WH_is_nondecreasing(input.recdate)) {
    // std::cout << "[WH_look_survdata] recdate = " << input.recdate.transpose() << std::endl;
    throw std::invalid_argument("[WH_look_survdata] recdates must be non-decreasing");
  }
  WH_survdata output;
  int n = (input.recdate.array() <= date).count();
  if (n == 0) {return output;}
  // Eigen::VectorXd output_recdate = input.recdate.head(n);
  // Eigen::VectorXi output_arm = input.arm.head(n);
  Eigen::VectorXd output_survdate = input.survdate.head(n);
  Eigen::VectorXd output_survstatus = input.survstatus.head(n);
  for (int i=0; i < n; i++) {
    if (output_survdate[i] > date) {
      output_survdate[i] = date;
      output_survstatus[i] = 0;
    }
  }
  output.recdate = input.recdate.head(n); // assumes sorted by recdate
  output.arm = input.arm.head(n);
  output.survdate = output_survdate;
  output.survtime = output_survdate - input.recdate.head(n);
  output.survstatus = output_survstatus;
  if (input.cluster.size() == input.recdate.size()) {
    output.cluster = input.cluster.head(n);
  }
  return output;
}

#endif
