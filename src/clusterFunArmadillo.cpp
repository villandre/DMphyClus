// [[Rcpp::depends(RcppArmadillo)]]

//#include <gperftools/profiler.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <algorithm>
#include "AugTree.h"
#include <limits>
#include <gsl/gsl_rng.h>


// [[Rcpp::plugins(openmp)]]

using namespace arma;
using namespace Rcpp;

template<class Mat>
void print_matrix(Mat matrix) {

    matrix.print(std::cout);
}

template<class Col>
void print_vector(Col colvec) {

    colvec.print(std::cout);
}

//provide explicit instantiations of the template function for
//every matrix type you use somewhere in your program.
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);
template void print_vector<arma::uvec>(arma::uvec colvec);
template void print_vector<arma::vec>(arma::vec colvec);

// [[Rcpp::export]]

List logLikCpp(IntegerMatrix & edgeMat, NumericVector & clusterMRCAs, NumericVector & limProbsVec, List & withinTransMatList, List & betweenTransMatList, int & numThreads, List & alignmentBin, uint & numTips, uint & numLoci, uint & withinMatListIndex, uint & betweenMatListIndex)
{
  //omp_set_num_threads(numOpenMP) ;
  std::vector<mat> withinTransMats = as<std::vector<mat>>(withinTransMatList) ;
  std::vector<mat> betweenTransMats = as<std::vector<mat>>(betweenTransMatList) ;
  vec limProbsVals = as<vec>(limProbsVec) ;
  
  std::vector<std::vector<uvec>> alignmentBinRecast = as<std::vector<std::vector<uvec>>>(alignmentBin) ;
  std::vector<std::vector<uvec>> * convertedBinData = new std::vector<std::vector<uvec>> ;
  convertedBinData->resize(alignmentBinRecast.size()) ;
  for (uint i = 0 ; i < alignmentBinRecast.size() ; i++) {
    convertedBinData->at(i).resize(alignmentBinRecast.at(i).size()) ;
    std::copy(alignmentBinRecast.at(i).begin(), alignmentBinRecast.at(i).end(), convertedBinData->at(i).begin()) ;
  }
  solutionDictionaryType solutionDictionary = new solutionDictionaryTypeNoPoint ;
  gsl_rng * randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  
  AugTree * PhylogeniesPoint1 = new AugTree(as<umat>(edgeMat), as<uvec>(clusterMRCAs), convertedBinData, solutionDictionary, withinMatListIndex, betweenMatListIndex, withinTransMatList.size(), randomNumGenerator, numThreads) ;

  PhylogeniesPoint1->ComputeLoglik(withinTransMats, betweenTransMats, limProbsVals) ;
  
  XPtr<AugTree> p(PhylogeniesPoint1, false) ; // Disabled automatic garbage collection. Tested with Valgrind, and no ensuing memory leak.
  XPtr<std::vector<std::vector<uvec>>> alignExtPoint(convertedBinData, false) ;
  
  return List::create(Named("logLik") = PhylogeniesPoint1->GetLoglik(),
                      Named("solutionPointer") = p, Named("alignmentBinPointer") = alignExtPoint) ;
}

std::unordered_map<std::string, uvec> defineMap(std::vector<std::string> & equivalency) 
{
  uint numStates = equivalency.size() ;
  std::unordered_map<std::string, uvec> myMap ;
  for (uint i = 0; i < numStates; i++) 
  {
    uvec unitVec(numStates, fill::zeros) ;
    unitVec(i) = 1 ;
    std::string label = equivalency.at(i) ;
    myMap[label] = unitVec ;
  }
  if (numStates == 4) 
  { // If we have four alphanumeric states, the assumption is that we're looking at a nucleotide alignment with states "a", "t", "c", "g".
    // Now, we handle ambiguities.
    myMap["r"] = (myMap["a"] + myMap["g"]) ;
    myMap["y"] = (myMap["c"] + myMap["t"]) ;
    myMap["s"] = (myMap["c"] + myMap["g"]) ;
    myMap["w"] = (myMap["a"] + myMap["t"]) ;
    myMap["k"] = (myMap["g"] + myMap["t"]) ;
    myMap["m"] = (myMap["a"] + myMap["c"]) ;
    myMap["b"] = (myMap["c"] + myMap["g"] + myMap["t"]) ;
    myMap["d"] = (myMap["a"] + myMap["g"] + myMap["t"]) ;
    myMap["h"] = (myMap["a"] + myMap["c"] + myMap["t"]) ;
    myMap["v"] = (myMap["a"] + myMap["c"] + myMap["g"]) ;
    myMap["-"] = (myMap["a"] + myMap["t"] + myMap["c"] + myMap["g"]) ;
    myMap["n"] = (myMap["a"] + myMap["t"] + myMap["c"] + myMap["g"]) ;
    myMap["."] = (myMap["a"] + myMap["t"] + myMap["c"] + myMap["g"]) ;
  }
  return myMap ;
}

// [[Rcpp::export]]

SEXP getConvertedAlignment(SEXP & equivVector, CharacterMatrix & alignmentAlphaMat)
{
  std::vector<std::vector<uvec>> alignmentBin ;
  alignmentBin.resize(alignmentAlphaMat.nrow()) ;
  std::vector<std::string> equivVectorAlpha = as<std::vector<std::string>>(equivVector) ;
  std::unordered_map<std::string, uvec> siteMap = defineMap(equivVectorAlpha) ;

  for (uint i = 0; i < alignmentAlphaMat.nrow(); i++) // Could probably be done more elegantly, i.e. with C++11 syntax...
  {
    alignmentBin.at(i).reserve(alignmentAlphaMat.ncol()) ;
    for (uint j = 0; j < alignmentAlphaMat.ncol(); j++) 
    {
      alignmentBin.at(i).push_back(siteMap[(std::string) alignmentAlphaMat(i,j)]) ;
    }
  }
  return wrap(alignmentBin) ;
}

// [[Rcpp::export]]

double newBetweenTransProbsLogLik(SEXP AugTreePointer, List & withinTransProbs, List & newBetweenTransProbs, IntegerMatrix & edgeMat, uint & numOpenMP, uint & newBetweenMatListIndex, NumericVector & limProbs) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(AugTreePointer == NULL)) 
  {
    XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
    std::vector<mat> withinTransMats = as<std::vector<mat>>(withinTransProbs) ;
    std::vector<mat> betweenTransMats = as<std::vector<mat>>(newBetweenTransProbs) ;
    vec limProbsVals = as<vec>(limProbs) ;
    
    pointedTree->NegateAllUpdateFlags() ;
    pointedTree->SetBetweenMatListIndex(newBetweenMatListIndex) ;
    pointedTree->CheckAndInvalidateBetweenRecursive(pointedTree->GetVertexVector().at(pointedTree->GetNumTips())) ;
    pointedTree->ComputeLoglik(withinTransMats, betweenTransMats, limProbsVals) ;
    
    return pointedTree->GetLoglik() ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

double newWithinTransProbsLogLik(SEXP AugTreePointer, List & newWithinTransProbs, List & betweenTransProbs, IntegerMatrix & edgeMat, uint & numOpenMP, uint & newWithinMatListIndex, NumericVector & limProbs) 
{
  //omp_set_num_threads(numOpenMP) ; 
  if (!(AugTreePointer == NULL)) 
  {
    XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
    std::vector<mat> withinTransMats = as<std::vector<mat>>(newWithinTransProbs) ;
    std::vector<mat> betweenTransMats = as<std::vector<mat>>(betweenTransProbs) ;
    vec limProbsVals = as<vec>(limProbs) ;
    
    pointedTree->NegateAllUpdateFlags() ;
    pointedTree->SetWithinMatListIndex(newWithinMatListIndex) ;
    
    pointedTree->InvalidateAll() ;
    pointedTree->ComputeLoglik(withinTransMats, betweenTransMats, limProbsVals) ;
    
    return pointedTree->GetLoglik() ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List withinClusNNIlogLik(SEXP AugTreePointer, IntegerMatrix & edgeMat, List & withinTransProbs, List & betweenTransProbs, uint & MRCAofClusForNNI, uint & numMovesNNI, uint & numOpenMP, NumericVector & limProbs) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(AugTreePointer == NULL)) 
  {
    std::vector<uint> vertexIndexForNNI ;
    std::vector<uint> vertexIndexVec ;
    
    XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
    std::vector<mat> withinTransMats = as<std::vector<mat>>(withinTransProbs) ;
    std::vector<mat> betweenTransMats = as<std::vector<mat>>(betweenTransProbs) ;
    vec limProbsVals = as<vec>(limProbs) ;
    
    pointedTree->NegateAllUpdateFlags() ;
    
    vertexIndexForNNI = pointedTree->GetNNIverticesWithin(pointedTree->GetVertexVector().at(MRCAofClusForNNI - 1)) ;
    
    for (uint counter = 0; counter < numMovesNNI; counter++)
    {
      unsigned long int rootForNNIindex = gsl_rng_uniform_int(pointedTree->GetRandomNumGenerator(), vertexIndexForNNI.size()) ;
      uvec placeholder(1) ;
      placeholder.at(0) = -1 ; // Probably a better way to do this than using a placeholder...
      vertexIndexVec = pointedTree->GetTwoVerticesForNNI(pointedTree->GetVertexVector().at(vertexIndexForNNI.at(rootForNNIindex)), placeholder) ;
      pointedTree->RearrangeTreeNNI(vertexIndexVec.at(0), vertexIndexVec.at(1)) ;
    }
    pointedTree->ComputeLoglik(withinTransMats, betweenTransMats, limProbsVals) ;
    umat newEdge = pointedTree->BuildEdgeMatrix() ; // All trees in the forest have the same hierarchy, hence the need to get the structure for only one of them.
    
    return List::create(Named("logLik") = pointedTree->GetLoglik(),
                        Named("edge") = newEdge) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List betweenClusNNIlogLik(SEXP AugTreePointer, List & withinTransProbs, List & betweenTransProbs, NumericVector & clusterMRCAs, IntegerMatrix & edgeMat, uint & numMovesNNI, uint & numOpenMP, NumericVector & limProbs) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(AugTreePointer == NULL)) 
  {
    std::vector<uint> vertexIndexForNNI ;
    std::vector<uint> vertexIndexVec ;
    uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
    
    XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
    std::vector<mat> withinTransMats = as<std::vector<mat>>(withinTransProbs) ;
    std::vector<mat> betweenTransMats = as<std::vector<mat>>(betweenTransProbs) ;
    vec limProbsVals = as<vec>(limProbs) ;
    
    pointedTree->NegateAllUpdateFlags() ;
    
    vertexIndexForNNI = pointedTree->GetNNIverticesBetween(pointedTree->GetVertexVector().at(pointedTree->GetNumTips()), clusterMRCAsRecast) ;
    
    for (uint counter = 0; counter < numMovesNNI; counter++)
    {
      unsigned long int rootForNNIindex = gsl_rng_uniform_int(pointedTree->GetRandomNumGenerator(), vertexIndexForNNI.size()) ;
      vertexIndexVec = pointedTree->GetTwoVerticesForNNI(pointedTree->GetVertexVector().at(vertexIndexForNNI.at(rootForNNIindex)), clusterMRCAsRecast) ;
      
      pointedTree->RearrangeTreeNNI(vertexIndexVec.at(0), vertexIndexVec.at(1)) ;
    }
    
    pointedTree->ComputeLoglik(withinTransMats, betweenTransMats, limProbsVals) ;
    umat newEdge = pointedTree->BuildEdgeMatrix() ;
    return List::create(Named("logLik") = pointedTree->GetLoglik(),
                       Named("edge") = newEdge) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

double clusSplitMergeLogLik(SEXP AugTreePointer, List & withinTransProbs, List & betweenTransProbs, IntegerVector & clusMRCAsToSplitOrMerge, IntegerMatrix & edgeMat, uint & numOpenMP, NumericVector & limProbs) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(AugTreePointer == NULL))
  {
    XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
    std::vector<mat> withinTransMats = as<std::vector<mat>>(withinTransProbs) ;
    std::vector<mat> betweenTransMats = as<std::vector<mat>>(betweenTransProbs) ;
    vec limProbsVals = as<vec>(limProbs) ;
    
    pointedTree->NegateAllUpdateFlags() ;
    
    uvec clusMRCAsToSplitOrMergeRecast = as<uvec>(clusMRCAsToSplitOrMerge) ;
   
    if (clusMRCAsToSplitOrMergeRecast.size() == 1) // Split move
    {
      pointedTree->HandleSplit(clusMRCAsToSplitOrMergeRecast.at(0)) ;
    }
    else
    {
      pointedTree->HandleMerge(clusMRCAsToSplitOrMergeRecast) ;
    }
    pointedTree->ComputeLoglik(withinTransMats, betweenTransMats, limProbsVals) ;
    
    return pointedTree->GetLoglik() ;
  }
  else
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

void RestorePreviousConfig(SEXP AugTreePointer, IntegerMatrix & edgeMat, int & withinMatListIndex, int & betweenMatListIndex, bool NNImove, IntegerVector & clusterMRCAs, bool splitMergeMove)
{
  XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
  pointedTree->RestorePreviousConfig(edgeMat, NNImove, withinMatListIndex, betweenMatListIndex, clusterMRCAs, splitMergeMove) ;
}

// [[Rcpp::export]]

void finalDeallocate(SEXP AugTreePointer) // We need to explicitly deallocate the random number generator.
{
  XPtr<AugTree> pointedTree(AugTreePointer) ; // Becomes a regular pointer again.
  gsl_rng_free(pointedTree->GetRandomNumGenerator()) ;
  delete pointedTree->GetSolutionDictionary() ;
  delete pointedTree->GetAlignmentBinReference() ;
  pointedTree->EndIOserviceAndJoinAll() ;
}