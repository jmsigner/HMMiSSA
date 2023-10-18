// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ll_Rcpp
double ll_Rcpp(arma::mat allprobs, arma::mat gamma, arma::rowvec delta, int nSteps);
RcppExport SEXP _HMMiSSA_ll_Rcpp(SEXP allprobsSEXP, SEXP gammaSEXP, SEXP deltaSEXP, SEXP nStepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type allprobs(allprobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type nSteps(nStepsSEXP);
    rcpp_result_gen = Rcpp::wrap(ll_Rcpp(allprobs, gamma, delta, nSteps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HMMiSSA_ll_Rcpp", (DL_FUNC) &_HMMiSSA_ll_Rcpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_HMMiSSA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
