// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// simple_integral
double simple_integral(arma::mat& lmat, arma::vec tseq);
RcppExport SEXP _TDAkit_simple_integral(SEXP lmatSEXP, SEXP tseqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lmat(lmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tseq(tseqSEXP);
    rcpp_result_gen = Rcpp::wrap(simple_integral(lmat, tseq));
    return rcpp_result_gen;
END_RCPP
}
// compute_slicemean
arma::mat compute_slicemean(arma::cube& dcube);
RcppExport SEXP _TDAkit_compute_slicemean(SEXP dcubeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type dcube(dcubeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_slicemean(dcube));
    return rcpp_result_gen;
END_RCPP
}
// compute_slicewsum
arma::mat compute_slicewsum(arma::cube& dcube, arma::vec weight);
RcppExport SEXP _TDAkit_compute_slicewsum(SEXP dcubeSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type dcube(dcubeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_slicewsum(dcube, weight));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _TDAkit_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _TDAkit_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _TDAkit_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _TDAkit_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TDAkit_simple_integral", (DL_FUNC) &_TDAkit_simple_integral, 2},
    {"_TDAkit_compute_slicemean", (DL_FUNC) &_TDAkit_compute_slicemean, 1},
    {"_TDAkit_compute_slicewsum", (DL_FUNC) &_TDAkit_compute_slicewsum, 2},
    {"_TDAkit_rcpparma_hello_world", (DL_FUNC) &_TDAkit_rcpparma_hello_world, 0},
    {"_TDAkit_rcpparma_outerproduct", (DL_FUNC) &_TDAkit_rcpparma_outerproduct, 1},
    {"_TDAkit_rcpparma_innerproduct", (DL_FUNC) &_TDAkit_rcpparma_innerproduct, 1},
    {"_TDAkit_rcpparma_bothproducts", (DL_FUNC) &_TDAkit_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_TDAkit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}