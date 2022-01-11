// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// discretizeCpp
Rcpp::List discretizeCpp(arma::mat a, arma::mat b, arma::rowvec w, arma::uvec mask, arma::umat ij, int Dx, int Dy);
RcppExport SEXP _MultiFit_discretizeCpp(SEXP aSEXP, SEXP bSEXP, SEXP wSEXP, SEXP maskSEXP, SEXP ijSEXP, SEXP DxSEXP, SEXP DySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type ij(ijSEXP);
    Rcpp::traits::input_parameter< int >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< int >::type Dy(DySEXP);
    rcpp_result_gen = Rcpp::wrap(discretizeCpp(a, b, w, mask, ij, Dx, Dy));
    return rcpp_result_gen;
END_RCPP
}
// single_Fisher_test
Rcpp::List single_Fisher_test(arma::rowvec t, bool correct, bool ret_all_probs);
RcppExport SEXP _MultiFit_single_Fisher_test(SEXP tSEXP, SEXP correctSEXP, SEXP ret_all_probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type t(tSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type ret_all_probs(ret_all_probsSEXP);
    rcpp_result_gen = Rcpp::wrap(single_Fisher_test(t, correct, ret_all_probs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MultiFit_discretizeCpp", (DL_FUNC) &_MultiFit_discretizeCpp, 7},
    {"_MultiFit_single_Fisher_test", (DL_FUNC) &_MultiFit_single_Fisher_test, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MultiFit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
