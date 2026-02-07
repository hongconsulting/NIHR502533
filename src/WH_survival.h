// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_SURVIVAL_H
#define WH_SURVIVAL_H

// misc tools for survival analysis

#include "WH_constexpr.h"
// #include "WH_IO.h"

double WH_h_to_S(const double hazard, const double time) {
  return std::exp(-hazard * time);
}

Eigen::VectorXd WH_h_to_S(const Eigen::VectorXd& hazard, const double time) {
  return (-hazard.array() * time).exp();
}

double WH_med_to_h(const double median) {
  return std::log(2.0) / median;
}

double WH_medcure_to_h(double med, double pi) {
  return (std::log(2.0) - std::log(1.0 - pi)) / med;
} // exp-distributed times among the uncured, median for the full mixture

double WH_med_to_S(const double median, const double time) {
  return std::exp(std::log(0.5) * time / median);
}

Eigen::VectorXd WH_med_to_logS(const double median, const Eigen::VectorXd& time) {
  return std::log(0.5) * time.array() / median;
}

Eigen::VectorXd WH_med_to_h(const Eigen::VectorXd& median) {
  return std::log(2.0) / median.array();
}

double WH_S_to_h(const double proportion, const double time) {
  return -std::log(proportion) / time;
}

Eigen::VectorXd WH_S_to_h(const Eigen::VectorXd& proportion, const double time) {
  return -proportion.array().log() / time;
}

double WH_S_to_med(const double proportion, const double time) {
  return -std::log(0.5) / (-std::log(proportion) / time);
}

#endif
