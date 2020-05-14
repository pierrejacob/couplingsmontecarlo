// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// ising_gibbs_sweep_
IntegerMatrix ising_gibbs_sweep_(IntegerMatrix state, NumericVector proba_beta);
RcppExport SEXP _couplingsmontecarlo_ising_gibbs_sweep_(SEXP stateSEXP, SEXP proba_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type state(stateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proba_beta(proba_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(ising_gibbs_sweep_(state, proba_beta));
    return rcpp_result_gen;
END_RCPP
}
// fast_rmvnorm_
NumericMatrix fast_rmvnorm_(int nsamples, const NumericVector& mean, const NumericMatrix& covariance);
RcppExport SEXP _couplingsmontecarlo_fast_rmvnorm_(SEXP nsamplesSEXP, SEXP meanSEXP, SEXP covarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type covariance(covarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_rmvnorm_(nsamples, mean, covariance));
    return rcpp_result_gen;
END_RCPP
}
// fast_rmvnorm_cholesky_
NumericMatrix fast_rmvnorm_cholesky_(int nsamples, const NumericVector& mean, const Eigen::MatrixXd& cholesky);
RcppExport SEXP _couplingsmontecarlo_fast_rmvnorm_cholesky_(SEXP nsamplesSEXP, SEXP meanSEXP, SEXP choleskySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type cholesky(choleskySEXP);
    rcpp_result_gen = Rcpp::wrap(fast_rmvnorm_cholesky_(nsamples, mean, cholesky));
    return rcpp_result_gen;
END_RCPP
}
// fast_dmvnorm_
NumericVector fast_dmvnorm_(const NumericMatrix& x, const NumericVector& mean, const NumericMatrix& covariance);
RcppExport SEXP _couplingsmontecarlo_fast_dmvnorm_(SEXP xSEXP, SEXP meanSEXP, SEXP covarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type covariance(covarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_dmvnorm_(x, mean, covariance));
    return rcpp_result_gen;
END_RCPP
}
// fast_dmvnorm_cholesky_inverse_
NumericVector fast_dmvnorm_cholesky_inverse_(const NumericMatrix& x, const NumericVector& mean, const Eigen::MatrixXd& cholesky_inverse);
RcppExport SEXP _couplingsmontecarlo_fast_dmvnorm_cholesky_inverse_(SEXP xSEXP, SEXP meanSEXP, SEXP cholesky_inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type cholesky_inverse(cholesky_inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_dmvnorm_cholesky_inverse_(x, mean, cholesky_inverse));
    return rcpp_result_gen;
END_RCPP
}
// rmvnorm_max_coupling_
Rcpp::List rmvnorm_max_coupling_(const NumericVector& mu1, const NumericVector& mu2, const NumericMatrix& Sigma1, const NumericMatrix& Sigma2);
RcppExport SEXP _couplingsmontecarlo_rmvnorm_max_coupling_(SEXP mu1SEXP, SEXP mu2SEXP, SEXP Sigma1SEXP, SEXP Sigma2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu2(mu2SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Sigma1(Sigma1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Sigma2(Sigma2SEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm_max_coupling_(mu1, mu2, Sigma1, Sigma2));
    return rcpp_result_gen;
END_RCPP
}
// rmvnorm_max_coupling_cholesky
Rcpp::List rmvnorm_max_coupling_cholesky(const NumericVector& mu1, const NumericVector& mu2, const Eigen::MatrixXd& Cholesky1, const Eigen::MatrixXd& Cholesky2, const Eigen::MatrixXd& Cholesky_inverse1, const Eigen::MatrixXd& Cholesky_inverse2);
RcppExport SEXP _couplingsmontecarlo_rmvnorm_max_coupling_cholesky(SEXP mu1SEXP, SEXP mu2SEXP, SEXP Cholesky1SEXP, SEXP Cholesky2SEXP, SEXP Cholesky_inverse1SEXP, SEXP Cholesky_inverse2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu2(mu2SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Cholesky1(Cholesky1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Cholesky2(Cholesky2SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Cholesky_inverse1(Cholesky_inverse1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Cholesky_inverse2(Cholesky_inverse2SEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm_max_coupling_cholesky(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2));
    return rcpp_result_gen;
END_RCPP
}
// rmvnorm_reflection_max_coupling_
Rcpp::List rmvnorm_reflection_max_coupling_(const Eigen::VectorXd& mu1, const Eigen::VectorXd& mu2, const Eigen::MatrixXd& Sigma_chol, const Eigen::MatrixXd& inv_Sigma_chol);
RcppExport SEXP _couplingsmontecarlo_rmvnorm_reflection_max_coupling_(SEXP mu1SEXP, SEXP mu2SEXP, SEXP Sigma_cholSEXP, SEXP inv_Sigma_cholSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mu2(mu2SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sigma_chol(Sigma_cholSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type inv_Sigma_chol(inv_Sigma_cholSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm_reflection_max_coupling_(mu1, mu2, Sigma_chol, inv_Sigma_chol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_couplingsmontecarlo_ising_gibbs_sweep_", (DL_FUNC) &_couplingsmontecarlo_ising_gibbs_sweep_, 2},
    {"_couplingsmontecarlo_fast_rmvnorm_", (DL_FUNC) &_couplingsmontecarlo_fast_rmvnorm_, 3},
    {"_couplingsmontecarlo_fast_rmvnorm_cholesky_", (DL_FUNC) &_couplingsmontecarlo_fast_rmvnorm_cholesky_, 3},
    {"_couplingsmontecarlo_fast_dmvnorm_", (DL_FUNC) &_couplingsmontecarlo_fast_dmvnorm_, 3},
    {"_couplingsmontecarlo_fast_dmvnorm_cholesky_inverse_", (DL_FUNC) &_couplingsmontecarlo_fast_dmvnorm_cholesky_inverse_, 3},
    {"_couplingsmontecarlo_rmvnorm_max_coupling_", (DL_FUNC) &_couplingsmontecarlo_rmvnorm_max_coupling_, 4},
    {"_couplingsmontecarlo_rmvnorm_max_coupling_cholesky", (DL_FUNC) &_couplingsmontecarlo_rmvnorm_max_coupling_cholesky, 6},
    {"_couplingsmontecarlo_rmvnorm_reflection_max_coupling_", (DL_FUNC) &_couplingsmontecarlo_rmvnorm_reflection_max_coupling_, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_couplingsmontecarlo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
