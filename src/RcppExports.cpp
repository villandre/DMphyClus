// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLikCpp
List logLikCpp(IntegerMatrix& edgeMat, NumericVector& clusterMRCAs, NumericVector& limProbsVec, List& withinTransMatList, List& betweenTransMatList, unsigned int& numThreads, List& alignmentBin, uint& numTips, uint& numLoci, uint& withinMatListIndex, uint& betweenMatListIndex);
RcppExport SEXP DMphyClus_logLikCpp(SEXP edgeMatSEXP, SEXP clusterMRCAsSEXP, SEXP limProbsVecSEXP, SEXP withinTransMatListSEXP, SEXP betweenTransMatListSEXP, SEXP numThreadsSEXP, SEXP alignmentBinSEXP, SEXP numTipsSEXP, SEXP numLociSEXP, SEXP withinMatListIndexSEXP, SEXP betweenMatListIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type clusterMRCAs(clusterMRCAsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbsVec(limProbsVecSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransMatList(withinTransMatListSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransMatList(betweenTransMatListSEXP);
    Rcpp::traits::input_parameter< unsigned int& >::type numThreads(numThreadsSEXP);
    Rcpp::traits::input_parameter< List& >::type alignmentBin(alignmentBinSEXP);
    Rcpp::traits::input_parameter< uint& >::type numTips(numTipsSEXP);
    Rcpp::traits::input_parameter< uint& >::type numLoci(numLociSEXP);
    Rcpp::traits::input_parameter< uint& >::type withinMatListIndex(withinMatListIndexSEXP);
    Rcpp::traits::input_parameter< uint& >::type betweenMatListIndex(betweenMatListIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikCpp(edgeMat, clusterMRCAs, limProbsVec, withinTransMatList, betweenTransMatList, numThreads, alignmentBin, numTips, numLoci, withinMatListIndex, betweenMatListIndex));
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
// newBetweenTransProbsLogLik
double newBetweenTransProbsLogLik(SEXP AugTreePointer, List& withinTransProbs, List& newBetweenTransProbs, IntegerMatrix& edgeMat, uint& numOpenMP, uint& newBetweenMatListIndex, NumericVector& limProbs);
RcppExport SEXP DMphyClus_newBetweenTransProbsLogLik(SEXP AugTreePointerSEXP, SEXP withinTransProbsSEXP, SEXP newBetweenTransProbsSEXP, SEXP edgeMatSEXP, SEXP numOpenMPSEXP, SEXP newBetweenMatListIndexSEXP, SEXP limProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransProbs(withinTransProbsSEXP);
    Rcpp::traits::input_parameter< List& >::type newBetweenTransProbs(newBetweenTransProbsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< uint& >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< uint& >::type newBetweenMatListIndex(newBetweenMatListIndexSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbs(limProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(newBetweenTransProbsLogLik(AugTreePointer, withinTransProbs, newBetweenTransProbs, edgeMat, numOpenMP, newBetweenMatListIndex, limProbs));
    return rcpp_result_gen;
END_RCPP
}
// newWithinTransProbsLogLik
double newWithinTransProbsLogLik(SEXP AugTreePointer, List& newWithinTransProbs, List& betweenTransProbs, IntegerMatrix& edgeMat, uint& numOpenMP, uint& newWithinMatListIndex, NumericVector& limProbs);
RcppExport SEXP DMphyClus_newWithinTransProbsLogLik(SEXP AugTreePointerSEXP, SEXP newWithinTransProbsSEXP, SEXP betweenTransProbsSEXP, SEXP edgeMatSEXP, SEXP numOpenMPSEXP, SEXP newWithinMatListIndexSEXP, SEXP limProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    Rcpp::traits::input_parameter< List& >::type newWithinTransProbs(newWithinTransProbsSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransProbs(betweenTransProbsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< uint& >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< uint& >::type newWithinMatListIndex(newWithinMatListIndexSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbs(limProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(newWithinTransProbsLogLik(AugTreePointer, newWithinTransProbs, betweenTransProbs, edgeMat, numOpenMP, newWithinMatListIndex, limProbs));
    return rcpp_result_gen;
END_RCPP
}
// withinClusNNIlogLik
List withinClusNNIlogLik(SEXP AugTreePointer, IntegerMatrix& edgeMat, List& withinTransProbs, List& betweenTransProbs, uint& MRCAofClusForNNI, uint& numMovesNNI, uint& numOpenMP, NumericVector& limProbs);
RcppExport SEXP DMphyClus_withinClusNNIlogLik(SEXP AugTreePointerSEXP, SEXP edgeMatSEXP, SEXP withinTransProbsSEXP, SEXP betweenTransProbsSEXP, SEXP MRCAofClusForNNISEXP, SEXP numMovesNNISEXP, SEXP numOpenMPSEXP, SEXP limProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransProbs(withinTransProbsSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransProbs(betweenTransProbsSEXP);
    Rcpp::traits::input_parameter< uint& >::type MRCAofClusForNNI(MRCAofClusForNNISEXP);
    Rcpp::traits::input_parameter< uint& >::type numMovesNNI(numMovesNNISEXP);
    Rcpp::traits::input_parameter< uint& >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbs(limProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(withinClusNNIlogLik(AugTreePointer, edgeMat, withinTransProbs, betweenTransProbs, MRCAofClusForNNI, numMovesNNI, numOpenMP, limProbs));
    return rcpp_result_gen;
END_RCPP
}
// betweenClusNNIlogLik
List betweenClusNNIlogLik(SEXP AugTreePointer, List& withinTransProbs, List& betweenTransProbs, NumericVector& clusterMRCAs, IntegerMatrix& edgeMat, uint& numMovesNNI, uint& numOpenMP, NumericVector& limProbs);
RcppExport SEXP DMphyClus_betweenClusNNIlogLik(SEXP AugTreePointerSEXP, SEXP withinTransProbsSEXP, SEXP betweenTransProbsSEXP, SEXP clusterMRCAsSEXP, SEXP edgeMatSEXP, SEXP numMovesNNISEXP, SEXP numOpenMPSEXP, SEXP limProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransProbs(withinTransProbsSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransProbs(betweenTransProbsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type clusterMRCAs(clusterMRCAsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< uint& >::type numMovesNNI(numMovesNNISEXP);
    Rcpp::traits::input_parameter< uint& >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbs(limProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(betweenClusNNIlogLik(AugTreePointer, withinTransProbs, betweenTransProbs, clusterMRCAs, edgeMat, numMovesNNI, numOpenMP, limProbs));
    return rcpp_result_gen;
END_RCPP
}
// clusSplitMergeLogLik
double clusSplitMergeLogLik(SEXP AugTreePointer, List& withinTransProbs, List& betweenTransProbs, IntegerVector& clusMRCAsToSplitOrMerge, IntegerMatrix& edgeMat, uint& numOpenMP, NumericVector& limProbs);
RcppExport SEXP DMphyClus_clusSplitMergeLogLik(SEXP AugTreePointerSEXP, SEXP withinTransProbsSEXP, SEXP betweenTransProbsSEXP, SEXP clusMRCAsToSplitOrMergeSEXP, SEXP edgeMatSEXP, SEXP numOpenMPSEXP, SEXP limProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransProbs(withinTransProbsSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransProbs(betweenTransProbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type clusMRCAsToSplitOrMerge(clusMRCAsToSplitOrMergeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< uint& >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type limProbs(limProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(clusSplitMergeLogLik(AugTreePointer, withinTransProbs, betweenTransProbs, clusMRCAsToSplitOrMerge, edgeMat, numOpenMP, limProbs));
    return rcpp_result_gen;
END_RCPP
}
// RestorePreviousConfig
void RestorePreviousConfig(SEXP AugTreePointer, IntegerMatrix& edgeMat, int& withinMatListIndex, int& betweenMatListIndex, bool NNImove, IntegerVector& clusterMRCAs, bool splitMergeMove);
RcppExport SEXP DMphyClus_RestorePreviousConfig(SEXP AugTreePointerSEXP, SEXP edgeMatSEXP, SEXP withinMatListIndexSEXP, SEXP betweenMatListIndexSEXP, SEXP NNImoveSEXP, SEXP clusterMRCAsSEXP, SEXP splitMergeMoveSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< int& >::type withinMatListIndex(withinMatListIndexSEXP);
    Rcpp::traits::input_parameter< int& >::type betweenMatListIndex(betweenMatListIndexSEXP);
    Rcpp::traits::input_parameter< bool >::type NNImove(NNImoveSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type clusterMRCAs(clusterMRCAsSEXP);
    Rcpp::traits::input_parameter< bool >::type splitMergeMove(splitMergeMoveSEXP);
    RestorePreviousConfig(AugTreePointer, edgeMat, withinMatListIndex, betweenMatListIndex, NNImove, clusterMRCAs, splitMergeMove);
    return R_NilValue;
END_RCPP
}
// finalDeallocate
void finalDeallocate(SEXP AugTreePointer);
RcppExport SEXP DMphyClus_finalDeallocate(SEXP AugTreePointerSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AugTreePointer(AugTreePointerSEXP);
    finalDeallocate(AugTreePointer);
    return R_NilValue;
END_RCPP
}
