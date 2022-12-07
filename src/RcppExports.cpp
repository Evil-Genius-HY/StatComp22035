// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// PageRank
NumericVector PageRank(NumericMatrix MA, double d, int MaxIter, NumericVector IniVec);
RcppExport SEXP _StatComp22035_PageRank(SEXP MASEXP, SEXP dSEXP, SEXP MaxIterSEXP, SEXP IniVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MA(MASEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type MaxIter(MaxIterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type IniVec(IniVecSEXP);
    rcpp_result_gen = Rcpp::wrap(PageRank(MA, d, MaxIter, IniVec));
    return rcpp_result_gen;
END_RCPP
}
// gibbsC
NumericMatrix gibbsC(NumericVector mu, NumericVector sigma, double rho, NumericVector initial, int N);
RcppExport SEXP _StatComp22035_gibbsC(SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP initialSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsC(mu, sigma, rho, initial, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp22035_PageRank", (DL_FUNC) &_StatComp22035_PageRank, 4},
    {"_StatComp22035_gibbsC", (DL_FUNC) &_StatComp22035_gibbsC, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp22035(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
