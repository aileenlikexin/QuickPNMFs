// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// PNMF_EucDistC
Eigen::MatrixXd PNMF_EucDistC(Eigen::MatrixXd X, Eigen::MatrixXd W_init, double tol, int maxIter, bool verboseN, double zerotol);
RcppExport SEXP _QuickPNMFs_PNMF_EucDistC(SEXP XSEXP, SEXP W_initSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP verboseNSEXP, SEXP zerotolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseN(verboseNSEXP);
    Rcpp::traits::input_parameter< double >::type zerotol(zerotolSEXP);
    rcpp_result_gen = Rcpp::wrap(PNMF_EucDistC(X, W_init, tol, maxIter, verboseN, zerotol));
    return rcpp_result_gen;
END_RCPP
}
// PNMF_KLC
Eigen::MatrixXd PNMF_KLC(Eigen::MatrixXd X, Eigen::MatrixXd W_init, double tol, int maxIter, bool verboseN, double zerotol);
RcppExport SEXP _QuickPNMFs_PNMF_KLC(SEXP XSEXP, SEXP W_initSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP verboseNSEXP, SEXP zerotolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseN(verboseNSEXP);
    Rcpp::traits::input_parameter< double >::type zerotol(zerotolSEXP);
    rcpp_result_gen = Rcpp::wrap(PNMF_KLC(X, W_init, tol, maxIter, verboseN, zerotol));
    return rcpp_result_gen;
END_RCPP
}
// DPNMFC
Eigen::MatrixXd DPNMFC(Eigen::MatrixXd X, Eigen::MatrixXd W_init, double tol, int maxIter, bool verboseN, double zerotol, Eigen::MatrixXd Xord, Eigen::VectorXi clunum, double mu, double lambda);
RcppExport SEXP _QuickPNMFs_DPNMFC(SEXP XSEXP, SEXP W_initSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP verboseNSEXP, SEXP zerotolSEXP, SEXP XordSEXP, SEXP clunumSEXP, SEXP muSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseN(verboseNSEXP);
    Rcpp::traits::input_parameter< double >::type zerotol(zerotolSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Xord(XordSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type clunum(clunumSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(DPNMFC(X, W_init, tol, maxIter, verboseN, zerotol, Xord, clunum, mu, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_QuickPNMFs_PNMF_EucDistC", (DL_FUNC) &_QuickPNMFs_PNMF_EucDistC, 6},
    {"_QuickPNMFs_PNMF_KLC", (DL_FUNC) &_QuickPNMFs_PNMF_KLC, 6},
    {"_QuickPNMFs_DPNMFC", (DL_FUNC) &_QuickPNMFs_DPNMFC, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_QuickPNMFs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
