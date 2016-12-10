// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLikCppToWrap
SEXP logLikCppToWrap(Rcpp::NumericMatrix& edgeMat, Rcpp::NumericMatrix& alignmentMat, Rcpp::NumericVector& limProbsVec, Rcpp::List& transMatList, int numOpenMP, SEXP& equivVector, Rcpp::CharacterMatrix& alignmentAlphaMat, Rcpp::List alignmentBin, Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag);
RcppExport SEXP DMphyClus_logLikCppToWrap(SEXP edgeMatSEXP, SEXP alignmentMatSEXP, SEXP limProbsVecSEXP, SEXP transMatListSEXP, SEXP numOpenMPSEXP, SEXP equivVectorSEXP, SEXP alignmentAlphaMatSEXP, SEXP alignmentBinSEXP, SEXP childNodeInClusIndicSEXP, SEXP returnMatIndicSEXP, SEXP internalFlagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type alignmentMat(alignmentMatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type limProbsVec(limProbsVecSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type transMatList(transMatListSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< SEXP& >::type equivVector(equivVectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterMatrix& >::type alignmentAlphaMat(alignmentAlphaMatSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type alignmentBin(alignmentBinSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type childNodeInClusIndic(childNodeInClusIndicSEXP);
    Rcpp::traits::input_parameter< const bool >::type returnMatIndic(returnMatIndicSEXP);
    Rcpp::traits::input_parameter< const bool >::type internalFlag(internalFlagSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikCppToWrap(edgeMat, alignmentMat, limProbsVec, transMatList, numOpenMP, equivVector, alignmentAlphaMat, alignmentBin, childNodeInClusIndic, returnMatIndic, internalFlag));
    return rcpp_result_gen;
END_RCPP
}
// getConvertedAlignmentToWrap
SEXP getConvertedAlignmentToWrap(Rcpp::NumericMatrix& alignmentMat, uint numStatesCons, int numOpenMP, SEXP& equivVector, Rcpp::CharacterMatrix& alignmentAlphaMat);
RcppExport SEXP DMphyClus_getConvertedAlignmentToWrap(SEXP alignmentMatSEXP, SEXP numStatesConsSEXP, SEXP numOpenMPSEXP, SEXP equivVectorSEXP, SEXP alignmentAlphaMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type alignmentMat(alignmentMatSEXP);
    Rcpp::traits::input_parameter< uint >::type numStatesCons(numStatesConsSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< SEXP& >::type equivVector(equivVectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterMatrix& >::type alignmentAlphaMat(alignmentAlphaMatSEXP);
    rcpp_result_gen = Rcpp::wrap(getConvertedAlignmentToWrap(alignmentMat, numStatesCons, numOpenMP, equivVector, alignmentAlphaMat));
    return rcpp_result_gen;
END_RCPP
}
// redimMultiBinByClus
SEXP redimMultiBinByClus(Rcpp::List multiBinByClus);
RcppExport SEXP DMphyClus_redimMultiBinByClus(SEXP multiBinByClusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type multiBinByClus(multiBinByClusSEXP);
    rcpp_result_gen = Rcpp::wrap(redimMultiBinByClus(multiBinByClus));
    return rcpp_result_gen;
END_RCPP
}
// logLikCppToWrapV
SEXP logLikCppToWrapV(Rcpp::List& edgeMatList, Rcpp::List& alignmentMatList, Rcpp::NumericVector& limProbsVec, Rcpp::List& transMatList, int numOpenMP, SEXP& equivVector, Rcpp::List& alignmentAlphaMatList, Rcpp::List alignmentBinList, Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const bool priorBySizeTransMatBool);
RcppExport SEXP DMphyClus_logLikCppToWrapV(SEXP edgeMatListSEXP, SEXP alignmentMatListSEXP, SEXP limProbsVecSEXP, SEXP transMatListSEXP, SEXP numOpenMPSEXP, SEXP equivVectorSEXP, SEXP alignmentAlphaMatListSEXP, SEXP alignmentBinListSEXP, SEXP childNodeInClusIndicSEXP, SEXP returnMatIndicSEXP, SEXP internalFlagSEXP, SEXP priorBySizeTransMatBoolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type edgeMatList(edgeMatListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type alignmentMatList(alignmentMatListSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type limProbsVec(limProbsVecSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type transMatList(transMatListSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< SEXP& >::type equivVector(equivVectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type alignmentAlphaMatList(alignmentAlphaMatListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type alignmentBinList(alignmentBinListSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type childNodeInClusIndic(childNodeInClusIndicSEXP);
    Rcpp::traits::input_parameter< const bool >::type returnMatIndic(returnMatIndicSEXP);
    Rcpp::traits::input_parameter< const bool >::type internalFlag(internalFlagSEXP);
    Rcpp::traits::input_parameter< const bool >::type priorBySizeTransMatBool(priorBySizeTransMatBoolSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikCppToWrapV(edgeMatList, alignmentMatList, limProbsVec, transMatList, numOpenMP, equivVector, alignmentAlphaMatList, alignmentBinList, childNodeInClusIndic, returnMatIndic, internalFlag, priorBySizeTransMatBool));
    return rcpp_result_gen;
END_RCPP
}
