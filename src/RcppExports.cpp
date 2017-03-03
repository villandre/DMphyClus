// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLikCpp
List logLikCpp(IntegerMatrix& edgeMat, NumericVector& clusterMRCAs, NumericVector& limProbsVec, List& withinTransMatList, List& betweenTransMatList, int numOpenMP, List alignmentBin, uint numTips, uint numLoci, uint withinMatListIndex, uint betweenMatListIndex);
RcppExport SEXP DMphyClus_logLikCpp(SEXP edgeMatSEXP, SEXP clusterMRCAsSEXP, SEXP limProbsVecSEXP, SEXP withinTransMatListSEXP, SEXP betweenTransMatListSEXP, SEXP numOpenMPSEXP, SEXP alignmentBinSEXP, SEXP numTipsSEXP, SEXP numLociSEXP, SEXP withinMatListIndexSEXP, SEXP betweenMatListIndexSEXP) {
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
    Rcpp::traits::input_parameter< uint >::type numTips(numTipsSEXP);
    Rcpp::traits::input_parameter< uint >::type numLoci(numLociSEXP);
    Rcpp::traits::input_parameter< uint >::type withinMatListIndex(withinMatListIndexSEXP);
    Rcpp::traits::input_parameter< uint >::type betweenMatListIndex(betweenMatListIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikCpp(edgeMat, clusterMRCAs, limProbsVec, withinTransMatList, betweenTransMatList, numOpenMP, alignmentBin, numTips, numLoci, withinMatListIndex, betweenMatListIndex));
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
// getNULLextPointer
SEXP getNULLextPointer();
RcppExport SEXP DMphyClus_getNULLextPointer() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getNULLextPointer());
    return rcpp_result_gen;
END_RCPP
}
// newBetweenTransProbsLogLik
List newBetweenTransProbsLogLik(SEXP ForestPointer, List& newBetweenTransProbs, IntegerMatrix& edgeMat, int numOpenMP, uint newBetweenMatListIndex);
RcppExport SEXP DMphyClus_newBetweenTransProbsLogLik(SEXP ForestPointerSEXP, SEXP newBetweenTransProbsSEXP, SEXP edgeMatSEXP, SEXP numOpenMPSEXP, SEXP newBetweenMatListIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    Rcpp::traits::input_parameter< List& >::type newBetweenTransProbs(newBetweenTransProbsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< uint >::type newBetweenMatListIndex(newBetweenMatListIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(newBetweenTransProbsLogLik(ForestPointer, newBetweenTransProbs, edgeMat, numOpenMP, newBetweenMatListIndex));
    return rcpp_result_gen;
END_RCPP
}
// newWithinTransProbsLogLik
List newWithinTransProbsLogLik(SEXP ForestPointer, List newWithinTransProbs, IntegerVector clusterMRCAs, IntegerMatrix& edgeMat, int numOpenMP, uint newWithinMatListIndex);
RcppExport SEXP DMphyClus_newWithinTransProbsLogLik(SEXP ForestPointerSEXP, SEXP newWithinTransProbsSEXP, SEXP clusterMRCAsSEXP, SEXP edgeMatSEXP, SEXP numOpenMPSEXP, SEXP newWithinMatListIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    Rcpp::traits::input_parameter< List >::type newWithinTransProbs(newWithinTransProbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type clusterMRCAs(clusterMRCAsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    Rcpp::traits::input_parameter< uint >::type newWithinMatListIndex(newWithinMatListIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(newWithinTransProbsLogLik(ForestPointer, newWithinTransProbs, clusterMRCAs, edgeMat, numOpenMP, newWithinMatListIndex));
    return rcpp_result_gen;
END_RCPP
}
// withinClusNNIlogLik
List withinClusNNIlogLik(SEXP ForestPointer, IntegerMatrix& edgeMat, uint MRCAofClusForNNI, uint numMovesNNI, int numOpenMP);
RcppExport SEXP DMphyClus_withinClusNNIlogLik(SEXP ForestPointerSEXP, SEXP edgeMatSEXP, SEXP MRCAofClusForNNISEXP, SEXP numMovesNNISEXP, SEXP numOpenMPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< uint >::type MRCAofClusForNNI(MRCAofClusForNNISEXP);
    Rcpp::traits::input_parameter< uint >::type numMovesNNI(numMovesNNISEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    rcpp_result_gen = Rcpp::wrap(withinClusNNIlogLik(ForestPointer, edgeMat, MRCAofClusForNNI, numMovesNNI, numOpenMP));
    return rcpp_result_gen;
END_RCPP
}
// betweenClusNNIlogLik
List betweenClusNNIlogLik(SEXP ForestPointer, NumericVector& clusterMRCAs, IntegerMatrix& edgeMat, uint numMovesNNI, int numOpenMP);
RcppExport SEXP DMphyClus_betweenClusNNIlogLik(SEXP ForestPointerSEXP, SEXP clusterMRCAsSEXP, SEXP edgeMatSEXP, SEXP numMovesNNISEXP, SEXP numOpenMPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type clusterMRCAs(clusterMRCAsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< uint >::type numMovesNNI(numMovesNNISEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    rcpp_result_gen = Rcpp::wrap(betweenClusNNIlogLik(ForestPointer, clusterMRCAs, edgeMat, numMovesNNI, numOpenMP));
    return rcpp_result_gen;
END_RCPP
}
// clusSplitMergeLogLik
List clusSplitMergeLogLik(SEXP ForestPointer, IntegerVector& clusMRCAsToSplitOrMerge, List& withinTransProbsMats, List& betweenTransProbsMats, IntegerMatrix& edgeMat, int numOpenMP);
RcppExport SEXP DMphyClus_clusSplitMergeLogLik(SEXP ForestPointerSEXP, SEXP clusMRCAsToSplitOrMergeSEXP, SEXP withinTransProbsMatsSEXP, SEXP betweenTransProbsMatsSEXP, SEXP edgeMatSEXP, SEXP numOpenMPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type clusMRCAsToSplitOrMerge(clusMRCAsToSplitOrMergeSEXP);
    Rcpp::traits::input_parameter< List& >::type withinTransProbsMats(withinTransProbsMatsSEXP);
    Rcpp::traits::input_parameter< List& >::type betweenTransProbsMats(betweenTransProbsMatsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type edgeMat(edgeMatSEXP);
    Rcpp::traits::input_parameter< int >::type numOpenMP(numOpenMPSEXP);
    rcpp_result_gen = Rcpp::wrap(clusSplitMergeLogLik(ForestPointer, clusMRCAsToSplitOrMerge, withinTransProbsMats, betweenTransProbsMats, edgeMat, numOpenMP));
    return rcpp_result_gen;
END_RCPP
}
// finalDeallocate
void finalDeallocate(SEXP ForestPointer);
RcppExport SEXP DMphyClus_finalDeallocate(SEXP ForestPointerSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    finalDeallocate(ForestPointer);
    return R_NilValue;
END_RCPP
}
// manualDeallocation
void manualDeallocation(SEXP ForestPointer);
RcppExport SEXP DMphyClus_manualDeallocation(SEXP ForestPointerSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ForestPointer(ForestPointerSEXP);
    manualDeallocation(ForestPointer);
    return R_NilValue;
END_RCPP
}
