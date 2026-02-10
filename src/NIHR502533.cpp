#include <Rcpp.h>

#include <string>
#include "WH_analyse_fixedt_summary.h"
#include "WH_binom_test.h"
#include "WH_concat.h"
#include "WH_constexpr.h"
#include "WH_distr_Bernoulli.h"
#include "WH_gen_survdataset.h"
// #include "WH_logrank_1sample.h"
#include "WH_look_survdata.h"
#include "WH_replicate.h"
#include "WH_string_f.h"
#include "WH_survdata.h"
// #include "WH_survdata_Rcpp.h"
#include "WH_survfit.h"
// #include "WH_survfit_Rcpp.h"
#include "WH_survival.h"

double Rcpp_BuyseTest(const Eigen::MatrixXd& outcomes, const Eigen::VectorXd& arm) {
  static bool loaded = false;
  if (!loaded) {
    Rcpp::Function sink("sink");
    Rcpp::Function library("library");
    sink(R_NilValue, Rcpp::Named("type") = "message");
    library("BuyseTest");
    sink(R_NilValue, Rcpp::Named("type") = "message");
    loaded = true;
  }
  Rcpp::Environment env = Rcpp::Environment::namespace_env("BuyseTest");
  Rcpp::Function BuyseTest = env["BuyseTest"];
  int n = outcomes.rows();
  int m = outcomes.cols();
  Rcpp::CharacterVector endpoint(m);
  Rcpp::List df(m + 1);
  df[0] = arm;
  for (int i = 0; i < m; ++i) {
    endpoint[i] = "V" + std::to_string(i + 1);
    df[i + 1] = outcomes.col(i);
  }
  Rcpp::CharacterVector colnames(m + 1);
  colnames[0] = "treatment";
  for (int i = 0; i < m; ++i) colnames[i + 1] = endpoint[i];
  df.attr("names") = colnames;
  df.attr("class") = "data.frame";
  df.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n);
  Rcpp::Function summary("summary");
  Rcpp::List s = summary(BuyseTest(
    Rcpp::Named("treatment") = "treatment",
    Rcpp::Named("data") = df,
    Rcpp::Named("endpoint") = endpoint,
    Rcpp::Named("type") = Rcpp::CharacterVector(m, "bin"),
    Rcpp::Named("operator") = Rcpp::CharacterVector(m, ">0"),
    Rcpp::Named("method.inference") = "varExact permutation",
    Rcpp::Named("trace") = 0
  ), Rcpp::Named("print") = false);
  Rcpp::String p = Rcpp::as<Rcpp::CharacterVector>(s["p.value"])[2];
  return std::stod(p.get_cstring());
}

Eigen::VectorXd WH_NIHR502533a_sim(const double n_target = 250,
                                      const double p_arm0surv = 0.73, const double p_arm1surv = 0.70, const double t_surv = 36,
                                      const double p_arm0tox1 = 0.24, const double p_arm1tox1 = 0.1,
                                      const double p_arm0tox2 = 0.45, const double p_arm1tox2 = 0.05,
                                      const double p_LTFU = 0.1, const double t_LTFU = 60,
                                      const double recrate = 7.5, const double t_recmax = 60, const double t_add = 36, const double t_max = 72) { // all arguments in months
  double t_add_d = t_add * WH::_DAYSPERMONTH; // all internal calculations in days
  double t_LTFU_d = t_LTFU * WH::_DAYSPERMONTH;
  double t_surv_d = t_surv * WH::_DAYSPERMONTH;
  double t_recmax_d = t_recmax * WH::_DAYSPERMONTH;
  double t_max_d = t_max * WH::_DAYSPERMONTH;

  Eigen::Vector2d lambda(WH_S_to_h(p_arm0surv, t_surv_d), WH_S_to_h(p_arm1surv, t_surv_d));
  Eigen::VectorXd lambda_LTFU = WH_replicate(WH_S_to_h(1.0 - p_LTFU, t_LTFU_d), 2);
  double lambda_rec = recrate / WH::_DAYSPERMONTH;
  Eigen::Vector2i rand_block(0, 1);
  Eigen::VectorXd rand_weight = Eigen::VectorXd::Zero(1);

  WH_survdata survdata = WH_gen_survdataset(lambda, lambda_rec, lambda_LTFU, n_target, rand_block, rand_weight, t_recmax_d, t_max_d);
  double t_newmax_d = std::min(survdata.recdate.maxCoeff() + t_add_d, t_max_d);

  WH_survdata look = WH_look_survdata(t_newmax_d, survdata);
  WH_survfit FT = WH_analyse_fixedt_summary(2, t_surv_d, look);
  double n_arm0survive = FT.risk[0];
  double n_arm0 = n_arm0survive + FT.event[0];
  double n_arm1survive = FT.risk[1];
  double n_arm1 = n_arm1survive + FT.event[1];
  int n = n_arm0 + n_arm1;

  Eigen::VectorXd arm0survive = Eigen::VectorXd::Zero(n_arm0); // for all outcomes, 1 is good
  arm0survive.head(n_arm0survive).setOnes();
  Eigen::VectorXd p_arm0notox1 = WH_replicate(1.0 - p_arm0tox1, n_arm0);
  Eigen::VectorXd arm0notox1 = WH_rBernoulli(p_arm0notox1);
  Eigen::VectorXd p_arm0notox2 = WH_replicate(1.0 - p_arm0tox2, n_arm0);
  Eigen::VectorXd arm0notox2 = WH_rBernoulli(p_arm0notox2);

  Eigen::VectorXd arm1survive = Eigen::VectorXd::Zero(n_arm1);
  arm1survive.head(n_arm1survive).setOnes();
  Eigen::VectorXd p_arm1notox1 = WH_replicate(1.0 - p_arm1tox1, n_arm1);
  Eigen::VectorXd arm1notox1 = WH_rBernoulli(p_arm1notox1);
  Eigen::VectorXd p_arm1notox2 = WH_replicate(1.0 - p_arm1tox2, n_arm1);
  Eigen::VectorXd arm1notox2 = WH_rBernoulli(p_arm1notox2);

  Eigen::VectorXd GPCarm = Eigen::VectorXd::Zero(n);
  GPCarm.tail(n_arm1).setOnes();

  Eigen::VectorXd GPCsurvive = WH_concat(arm0survive, arm1survive);
  Eigen::VectorXd GPCnotox1 = WH_concat(arm0notox1, arm1notox1);
  Eigen::VectorXd GPCnotox2 = WH_concat(arm0notox2, arm1notox2);

  Eigen::MatrixXd GPCoutcomes = Eigen::MatrixXd(n, 3);
  GPCoutcomes.col(0) = GPCsurvive;
  GPCoutcomes.col(1) = GPCnotox1;
  GPCoutcomes.col(2) = GPCnotox2;

  // Eigen::MatrixXd GPC = WH_GPC(GPCoutcomes, GPCarm);
  // Eigen::VectorXi survstatus = look.survstatus.cast<int>();
  // Eigen::MatrixXd arm = look.arm.cast<double>();
  // WH_regfit Cox = WH_Cox(arm, look.survtime, survstatus);

  // Eigen::VectorXd output(6);
  // output[1] = GPC(2, 5); // NTB
  // output[2] = GPC(2, 7); // NTBp
  // output[3] = Cox.coef[0];
  // output[4] = Cox.p[0];
  // output[5] = t_newmax_d/365.2425;
  Eigen::VectorXd output(3);
  output[0] = n;
  output[1] = Rcpp_BuyseTest(GPCoutcomes, GPCarm);
  output[2] = t_newmax_d/365.2425;
  return output;
}

//' NIHR502533 cohort A
//'
//' Sample size simulator for NIHR502533 cohort A.
//'
//' @param name Scenario label. Default = `"scenario 1"`.
//' @param n_target Target sample size. Default = `204`.
//' @param p_arm0surv Control arm recurrence-free survival probability at `t_surv`
//' months. Default = `0.73`.
//' @param p_arm1surv Experimental arm recurrence-free survival probability at
//' `t_surv` months. Default = `0.70`.
//' @param t_surv Time in months when survival is assessed. Default = `36`.
//' @param p_arm0tox1 Control arm probability of toxicity 1. Default = `0.24`.
//' @param p_arm1tox1 Experimental arm probability of toxicity 1. Default =
//' `0.10`.
//' @param p_arm0tox2 Control arm probability of toxicity 2. Default = `0.45`.
//' @param p_arm1tox2 Experimental arm probability of toxicity 2. Default =
//' `0.05`.
//' @param p_LTFU Average proportion lost to follow-up per `t_LTFU` months.
//' Default = `0.1`.
//' @param t_LTFU Duration in months over which `p_LTFU` are lost to follow-up.
//' Default = `60`.
//' @param recrate Average number of participants recruited per month. Default
//' = `7.5`.
//' @param t_recmax Maximum recruitment duration in months. Default = `60`.
//' @param t_add Months of additional follow-up after recruitment closes.
//' Default = `36`.
//' @param t_max Maximum trial duration in months. Default = `72`.
//' @param reps Number of simulation replicates. Default = `1e5`.
//' @param seed Random seed. Default = `24601`.
//' @param echo Logical value indicating whether to print scenario information.
//' Default = `TRUE`.
//' @return A numeric vector with 3 elements: mean sample size, power, and mean
//' trial duration in months.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd WH_NIHR502533a(std::string name = "scenario 1", double n_target = 204,
                              double p_arm0surv = 0.73, double p_arm1surv = 0.70,
                              double t_surv = 36,
                              double p_arm0tox1 = 0.25, double p_arm1tox1 = 0.1,
                              double p_arm0tox2 = 0.45, double p_arm1tox2 = 0.05,
                              double p_LTFU = 0.1, double t_LTFU = 60,
                              double recrate = 7.5, double t_recmax = 60,
                              double t_add = 36, double t_max = 72,
                              int reps = 1e5, int seed = 24601, bool echo = true) {
  Eigen::MatrixXd sims = Eigen::MatrixXd(reps, 3);
  // int n_thread = omp_get_max_threads();
  // omp_set_num_threads(n_thread);
  // int chunksize = reps / n_thread;
  if (echo) {
    Rcpp::Rcout << "NIHR502533 cohort A: " << name << "\n" << " - Target N = " << n_target << "\n";
    Rcpp::Rcout << " - RFS at " << t_surv/12 << " years = " << 100 * p_arm0surv << "% vs " << 100 * p_arm1surv << "%" << "\n";
    Rcpp::Rcout << " - Tox 1 = " << 100 * p_arm0tox1 << "% vs " << 100 * p_arm1tox1 << "%" << "\n";
    Rcpp::Rcout << " - Tox 2 = " << 100 * p_arm0tox2 << "% vs " << 100 * p_arm1tox2 << "%" << "\n";
    Rcpp::Rcout << " - Average recruitment = " << recrate << "/month" << "\n";
    Rcpp::Rcout << " - Loss to follow-up = " << 100 * p_LTFU << "% every " << t_LTFU/12 << " years" << "\n";
    Rcpp::Rcout << " - Recruitment duration = up to " << t_recmax/12 << " years" << "\n";
    Rcpp::Rcout << " - Max duration = " << t_add/12 << " years after the last recruitment, up to " << t_max/12 << " years" << "\n";
  }
  // std::cout << " *** Simulation of " << WH_string_SI(reps, 1) << " trials (chunked parallel: " << n_thread << " threads) ***" << "\n";
  // #pragma omp parallel num_threads(n_thread)
  // {
  // WH_RNG.seed(seed + omp_get_thread_num());
  WH_RNG.seed(seed);
  // int i_thread = omp_get_thread_num();
  // int i_start = i_thread * chunksize;
  int i_start = 0;
  // int i_end = (i_thread == n_thread - 1) ? reps : (i_start + chunksize); // last thread takes leftovers
  int i_end = reps;
  for (int i = i_start; i < i_end; i++) {
    sims.row(i) = WH_NIHR502533a_sim(n_target, p_arm0surv, p_arm1surv, t_surv, p_arm0tox1, p_arm1tox1, p_arm0tox2, p_arm1tox2, p_LTFU,
             t_LTFU, recrate, t_recmax, t_add, t_max);
  }
  // }

  // Eigen::VectorXd output(6);
  // output[0] = sims.col(0).mean(); // n
  // output[1] = sims.col(1).mean(); // NTB
  // output[2] = (sims.col(2).array() < 0.05).cast<double>().mean();
  // output[3] = std::exp(sims.col(3).mean());
  // output[4] = (sims.col(4).array() < 0.05).cast<double>().mean();
  // output[5] = (sims.col(5).array().mean());
  Eigen::VectorXd output(3);
  output[0] = sims.col(0).mean(); // n
  output[1] = (sims.col(1).array() < 0.05).cast<double>().mean(); // power
  output[2] = (sims.col(2).array().mean()); // duration
  if (echo) {
    Rcpp::Rcout << " - Average N = " << WH_string_f(output[0]) << "\n";
    Rcpp::Rcout << " - Power (NTB) = " << WH_string_f(100 * output[1]) << "%" << "\n";
    Rcpp::Rcout << " - Average duration = " << WH_string_f(output[2]) << " years" << "\n";
  }
  // std::cout << " - Average N = " << WH_string_f(output[0]) << "\n";
  // std::cout << " - Average duration = " << WH_string_f(output[5]) << " years" << "\n";
  // std::cout << " - Average effect (NTB) = " << WH_string_f(output[1]);
  // std::cout << ", power (NTB) = " << WH_string_f(100 * output[2]) << "%" << "\n";
  // std::cout << " - Average effect (RFS HR) = " << WH_string_f(output[3]);
  // std::cout << ", power (RFS HR) = " << WH_string_f(100 * output[4]) << "%" << "\n" << "\n";
  return output;
}

Eigen::VectorXd WH_NIHR502533b_sim(const double n_target = 19,
                                          const double p_nullclear = 0.1,
                                          const double p_arm0clear = 0.35,
                                          const double t_surv = 12,
                                          const double p_LTFU = 0.15, const double t_LTFU = 60,
                                          const double recrate = 7.5, const double t_recmax = 60,
                                          const double t_add = 36, const double t_max = 72) { // all arguments in months
  double p_nullsurv = 1 - p_nullclear;
  double p_arm0surv = 1 - p_arm0clear;
  double t_add_d = t_add * WH::_DAYSPERMONTH; // all internal calculations in days
  double t_LTFU_d = t_LTFU * WH::_DAYSPERMONTH;
  double t_surv_d = t_surv * WH::_DAYSPERMONTH;
  double t_recmax_d = t_recmax * WH::_DAYSPERMONTH;
  double t_max_d = t_max * WH::_DAYSPERMONTH;

  Eigen::VectorXd lambda_null = Eigen::VectorXd::Constant(1, WH_S_to_h(p_nullsurv, t_surv_d));
  Eigen::VectorXd lambda = Eigen::VectorXd::Constant(1, WH_S_to_h(p_arm0surv, t_surv_d));
  Eigen::VectorXd lambda_LTFU = Eigen::VectorXd::Constant(1, WH_S_to_h(1.0 - p_LTFU, t_LTFU_d));
  double lambda_rec = recrate / WH::_DAYSPERMONTH;

  Eigen::VectorXi rand_block = Eigen::VectorXi::Zero(1);
  Eigen::VectorXd rand_weight = Eigen::VectorXd::Zero(1);

  WH_survdata survdata = WH_gen_survdataset(lambda, lambda_rec, lambda_LTFU, n_target, rand_block, rand_weight, t_recmax_d, t_max_d);
  double t_newmax_d = std::min(survdata.recdate.maxCoeff() + t_add_d, t_max_d);

  WH_survdata look = WH_look_survdata(t_newmax_d, survdata);
  WH_survfit FT = WH_analyse_fixedt_summary(1, t_surv_d, look);
  double n_arm0survive = FT.risk[0];
  double n_arm0 = n_arm0survive + FT.event[0];
  Eigen::VectorXd arm0survive = Eigen::VectorXd::Zero(n_arm0);
  arm0survive.head(n_arm0survive).setOnes();

  // Eigen::VectorXd output(5);
  // output[0] = n_arm0;
  // output[1] = n_arm0 - n_arm0survive; // number of clearances
  // output[2] = WH_logrank_1sample(look.survtime, look.survstatus, lambda_null[0]);
  // output[3] = WH_binom_test(n_arm0survive, n_arm0, p_nullsurv);
  // output[4] = t_newmax_d/365.2425;
  Eigen::VectorXd output(3);
  output[0] = n_arm0;
  output[1] = WH_binom_test(n_arm0survive, n_arm0, p_nullsurv);
  output[2] = t_newmax_d/365.2425;
  return output;
}

//' NIHR502533 cohort B
//'
//' Sample size simulator for NIHR502533 cohort B.
//'
//' @param name Scenario label. Default = `"scenario 1"`.
//' @param n_target Target sample size. Default = `19`.
//' @param p_nullclear Null hypothesis probability of clearance at `t_surv`
//' months. Default = `0.1`.
//' @param p_arm0clear Expected probability of clearance under the alternative
//' at `t_surv` months. Default = `0.35`.
//' @param t_surv Time in months when clearance is assessed. Default = `12`.
//' @param p_LTFU Average proportion lost to follow-up per `t_LTFU` months.
//' Default = `0.15`.
//' @param t_LTFU Duration in months over which `p_LTFU` are lost to follow-up.
//' Default = `60`.
//' @param recrate Average number of participants recruited per month. Default =
//' `7.5`.
//' @param t_recmax Maximum recruitment duration in months. Default = `60`.
//' @param t_add Months of additional follow-up after recruitment closes.
//' Default = `36`.
//' @param t_max Maximum trial duration in months. Default = `72`.
//' @param reps Number of simulation replicates. Default = `1e5`.
//' @param seed Random seed. Default = `24601`.
//' @return A numeric vector with 3 elements: mean sample size, power, and mean
//' trial duration in months.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd WH_NIHR502533b(std::string name = "scenario 1",
                                      double n_target = 19,
                                      const double p_nullclear = 0.1,
                                      const double p_arm0clear = 0.35,
                                      double t_surv = 12,
                                      const double p_LTFU = 0.15, const double t_LTFU = 60,
                                      const double recrate = 7.5, const double t_recmax = 60,
                                      double t_add = 36, double t_max = 72, int reps = 1e5, int seed = 24601) {
  Eigen::MatrixXd sims = Eigen::MatrixXd(reps, 3);
  // int n_thread = omp_get_max_threads();
  // omp_set_num_threads(n_thread);
  // int chunksize = reps / n_thread;
  Rcpp::Rcout << "NIHR502533 cohort B: " << name << std::endl << " - Target N = " << n_target << std::endl;
  Rcpp::Rcout << " - Clearance at " << t_surv/12 << " years = " << 100 * p_arm0clear << "%" << std::endl;
  Rcpp::Rcout << " - Average recruitment = " << recrate << "/month" << std::endl;
  Rcpp::Rcout << " - Loss to follow-up = " << 100 * p_LTFU << "% every " << t_LTFU/12 << " years" << std::endl;
  Rcpp::Rcout << " - Recruitment duration = up to " << t_recmax/12 << " years" << std::endl;
  Rcpp::Rcout << " - Max duration = " << t_add/12 << " years after the last recruitment, up to " << t_max/12 << " years" << std::endl;
  WH_RNG.seed(seed);
  // std::cout << " *** Simulation of " << WH_string_SI(reps, 1) << " trials (chunked parallel: " << n_thread << " threads) ***" << std::endl;
  // #pragma omp parallel num_threads(n_thread)
  // {
  // WH_RNG.seed(seed + omp_get_thread_num());
  //int i_thread = omp_get_thread_num();
  //int i_start = i_thread * chunksize;
  int i_start = 0;
  //int i_end = (i_thread == n_thread - 1) ? reps : (i_start + chunksize); // last thread takes leftovers
  int i_end = reps;
  for (int i = i_start; i < i_end; i++) {
    sims.row(i) = WH_NIHR502533b_sim(n_target, p_nullclear, p_arm0clear, t_surv,
             p_LTFU, t_LTFU, recrate, t_recmax, t_add, t_max);
  }
  // }
  // Eigen::VectorXd output(5);
  // output[0] = sims.col(0).mean(); // n
  // output[1] = sims.col(1).mean(); // n_cleared
  // output[2] = (sims.col(2).array() < 0.05).cast<double>().mean(); // logrank
  // output[3] = (sims.col(3).array() < 0.05).cast<double>().mean(); // binom
  // output[4] = (sims.col(4).array().mean()); // duration
  Eigen::VectorXd output(3);
  output[0] = sims.col(0).mean(); // n
  output[1] = (sims.col(1).array() < 0.05).cast<double>().mean(); // power
  output[2] = (sims.col(2).array().mean()); // duration
  Rcpp::Rcout << " - Average N = " << WH_string_f(output[0]) << std::endl;
  Rcpp::Rcout << " - Power = " << WH_string_f(100 * output[1]) << "%" << std::endl;
  Rcpp::Rcout << " - Average duration = " << WH_string_f(output[2]) << " years" << std::endl;
  return output;
}
