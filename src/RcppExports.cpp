// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// compute_log10BF
arma::colvec compute_log10BF(double bhat, double sdhat, arma::colvec& phi_vec);
RcppExport SEXP _BagseCpp_compute_log10BF(SEXP bhatSEXP, SEXP sdhatSEXP, SEXP phi_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type bhat(bhatSEXP);
    Rcpp::traits::input_parameter< double >::type sdhat(sdhatSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type phi_vec(phi_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_log10BF(bhat, sdhat, phi_vec));
    return rcpp_result_gen;
END_RCPP
}
// weightedlog10BF
double weightedlog10BF(arma::colvec& log10_BFv, arma::colvec& wv);
RcppExport SEXP _BagseCpp_weightedlog10BF(SEXP log10_BFvSEXP, SEXP wvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type log10_BFv(log10_BFvSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type wv(wvSEXP);
    rcpp_result_gen = Rcpp::wrap(weightedlog10BF(log10_BFv, wv));
    return rcpp_result_gen;
END_RCPP
}
// computegridwts
arma::colvec computegridwts(arma::colvec& log10_BF_vec, arma::colvec& pi_vec, double pi0, double log10_NC);
RcppExport SEXP _BagseCpp_computegridwts(SEXP log10_BF_vecSEXP, SEXP pi_vecSEXP, SEXP pi0SEXP, SEXP log10_NCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type log10_BF_vec(log10_BF_vecSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< double >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< double >::type log10_NC(log10_NCSEXP);
    rcpp_result_gen = Rcpp::wrap(computegridwts(log10_BF_vec, pi_vec, pi0, log10_NC));
    return rcpp_result_gen;
END_RCPP
}
// torus_pool_em
arma::colvec torus_pool_em(arma::colvec& p, arma::mat& BF_matrix, arma::uvec& annot_vec);
RcppExport SEXP _BagseCpp_torus_pool_em(SEXP pSEXP, SEXP BF_matrixSEXP, SEXP annot_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type BF_matrix(BF_matrixSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type annot_vec(annot_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(torus_pool_em(p, BF_matrix, annot_vec));
    return rcpp_result_gen;
END_RCPP
}
// torus_pool_loglik
double torus_pool_loglik(arma::vec& p, arma::mat& BF_matrix, arma::uvec& annot_vec);
RcppExport SEXP _BagseCpp_torus_pool_loglik(SEXP pSEXP, SEXP BF_matrixSEXP, SEXP annot_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type BF_matrix(BF_matrixSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type annot_vec(annot_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(torus_pool_loglik(p, BF_matrix, annot_vec));
    return rcpp_result_gen;
END_RCPP
}
// torus_cpp
List torus_cpp(arma::colvec& betahat, arma::colvec& sebetahat, arma::uvec& annotation, double tol);
RcppExport SEXP _BagseCpp_torus_cpp(SEXP betahatSEXP, SEXP sebetahatSEXP, SEXP annotationSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type sebetahat(sebetahatSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type annotation(annotationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(torus_cpp(betahat, sebetahat, annotation, tol));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _BagseCpp_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _BagseCpp_rcpparma_outerproduct(SEXP xSEXP) {
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
RcppExport SEXP _BagseCpp_rcpparma_innerproduct(SEXP xSEXP) {
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
RcppExport SEXP _BagseCpp_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BagseCpp_compute_log10BF", (DL_FUNC) &_BagseCpp_compute_log10BF, 3},
    {"_BagseCpp_weightedlog10BF", (DL_FUNC) &_BagseCpp_weightedlog10BF, 2},
    {"_BagseCpp_computegridwts", (DL_FUNC) &_BagseCpp_computegridwts, 4},
    {"_BagseCpp_torus_pool_em", (DL_FUNC) &_BagseCpp_torus_pool_em, 3},
    {"_BagseCpp_torus_pool_loglik", (DL_FUNC) &_BagseCpp_torus_pool_loglik, 3},
    {"_BagseCpp_torus_cpp", (DL_FUNC) &_BagseCpp_torus_cpp, 4},
    {"_BagseCpp_rcpparma_hello_world", (DL_FUNC) &_BagseCpp_rcpparma_hello_world, 0},
    {"_BagseCpp_rcpparma_outerproduct", (DL_FUNC) &_BagseCpp_rcpparma_outerproduct, 1},
    {"_BagseCpp_rcpparma_innerproduct", (DL_FUNC) &_BagseCpp_rcpparma_innerproduct, 1},
    {"_BagseCpp_rcpparma_bothproducts", (DL_FUNC) &_BagseCpp_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BagseCpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
