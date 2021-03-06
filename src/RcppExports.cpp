// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// center_cpp
SEXP center_cpp(Rcpp::NumericMatrix m, Rcpp::LogicalVector inplace);
RcppExport SEXP _ccf_center_cpp(SEXP mSEXP, SEXP inplaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type inplace(inplaceSEXP);
    rcpp_result_gen = Rcpp::wrap(center_cpp(m, inplace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ccf_center_cpp", (DL_FUNC) &_ccf_center_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ccf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
