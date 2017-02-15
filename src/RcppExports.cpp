// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLikCpp
List logLikCpp(IntegerMatrix& edgeMat, NumericVector& clusterMRCAs, NumericVector& limProbsVec, List& withinTransMatList, List& betweenTransMatList, int numOpenMP, List alignmentBin, int numTips);
RcppExport SEXP DMphyClus_logLikCpp(SEXP edgeMatSEXP, SEXP clusterMRCAsSEXP, SEXP limProbsVecSEXP, SEXP withinTransMatListSEXP, SEXP betweenTransMatListSEXP, SEXP numOpenMPSEXP, SEXP alignmentBinSEXP, SEXP numTipsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type clusterMRCAs(clusterMRCAsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbsVec(limProbsVecSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransMatList(withinTransMatListSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransMatList(betweenTransMatListSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< List >::type alignmentBin(alignmentBinSEXP);
    Rcpp::traits::input_parameter< int >::type numTips(numTipsSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikCpp(edgeMat, clusterMRCAs, limProbsVec, withinTransMatList, betweenTransMatList, numOpenMP, alignmentBin, numTips));
    return rcpp_result_gen;
END_RCPP
}
// getConvertedAlignment
SEXP getConvertedAlignment(SEXP& equivVector, CharacterMatrix& alignmentAlphaMat);
RcppExport SEXP DMphyClus_getConvertedAlignment(SEXP equivVectorSEXP, SEXP alignmentAlphaMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP& >::type equivVector(equivVectorSEXP);
    Rcpp::traits::input_parameter< CharacterMatrix& >::type alignmentAlphaMat(alignmentAlphaMatSEXP);
    rcpp_result_gen = Rcpp::wrap(getConvertedAlignment(equivVector, alignmentAlphaMat));
    return rcpp_result_gen;
END_RCPP
}
