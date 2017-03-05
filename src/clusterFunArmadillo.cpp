// [[Rcpp::depends(RcppArmadillo)]]

//#include <gperftools/profiler.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <algorithm>
#include "Forest.h"
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

List logLikCpp(IntegerMatrix & edgeMat, NumericVector & clusterMRCAs, NumericVector & limProbsVec, List & withinTransMatList, List & betweenTransMatList, int numOpenMP, List alignmentBin, uint numTips, uint numLoci, uint withinMatListIndex, uint betweenMatListIndex)
{
  //omp_set_num_threads(numOpenMP) ;
  
  std::vector<std::vector<uvec>> alignmentBinRecast = as<std::vector<std::vector<uvec>>>(alignmentBin) ;
  std::vector<std::vector<uvec>> * convertedBinData = new std::vector<std::vector<uvec>> ;
  convertedBinData->resize(alignmentBinRecast.size()) ;
  for (uint i = 0 ; i < alignmentBinRecast.size() ; i++) {
    convertedBinData->at(i).resize(alignmentBinRecast.at(i).size()) ;
    std::copy(alignmentBinRecast.at(i).begin(), alignmentBinRecast.at(i).end(), convertedBinData->at(i).begin()) ;
  }
  solutionDictionaryType solutionDictionary = new std::vector<std::unordered_map<std::size_t, vec>>(withinTransMatList.size()) ;
  Forest * PhylogeniesPoint1 = new Forest(edgeMat, clusterMRCAs, convertedBinData, withinTransMatList, betweenTransMatList, limProbsVec, numTips, numLoci, solutionDictionary, withinMatListIndex, betweenMatListIndex);
  PhylogeniesPoint1->ComputeLoglik() ;
  
  Forest * PhylogeniesPoint2 = new Forest(edgeMat, clusterMRCAs, convertedBinData, withinTransMatList, betweenTransMatList, limProbsVec, numTips, numLoci, solutionDictionary, withinMatListIndex, betweenMatListIndex);
  gsl_rng_free(PhylogeniesPoint2->GetRandomNumGenerator()) ; // We want the random number generator to be the same for both Forests, hence this deallocation, without which there would be a memory leak.
  PhylogeniesPoint2->SetRNG(PhylogeniesPoint1->GetRandomNumGenerator()) ; 
  
  XPtr<Forest> p(PhylogeniesPoint1, false) ; // Disabled automatic garbage collection. Tested with Valgrind, and no ensuing memory leak.
  XPtr<std::vector<std::vector<uvec>>> alignExtPoint(convertedBinData, false) ;
  XPtr<Forest> p2(PhylogeniesPoint2, false) ;
  
  return List::create(Named("logLik") = PhylogeniesPoint1->GetLoglik(),
                      Named("solutionPointer") = p, Named("alignmentBinPointer") = alignExtPoint, Named("alternatePointer") = p2) ;
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
  alignmentBin.resize(alignmentAlphaMat.ncol()) ;
  std::vector<std::string> equivVectorAlpha = as<std::vector<std::string>>(equivVector) ;
  std::unordered_map<std::string, uvec> siteMap = defineMap(equivVectorAlpha) ;

  for (uint i = 0; i < alignmentAlphaMat.ncol(); i++) // Could probably be done more elegantly, i.e. with C++11 syntax...
  {
    alignmentBin.at(i).reserve(alignmentAlphaMat.nrow()) ;
    for (uint j = 0; j < alignmentAlphaMat.nrow(); j++) 
    {
      alignmentBin.at(i).push_back(siteMap[(std::string) alignmentAlphaMat(j,i)]) ;
    }
  }
  return wrap(alignmentBin) ;
}

// [[Rcpp::export]]

SEXP getNULLextPointer()
{
  Rcpp::XPtr<int> p(NULL, true) ;
  return p ;
}

// [[Rcpp::export]]

List newBetweenTransProbsLogLik(SEXP ForestPointer, SEXP alternatePointer, List & newBetweenTransProbs, IntegerMatrix & edgeMat, uint & numOpenMP, uint & newBetweenMatListIndex) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointer == NULL)) 
  {
    XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
    XPtr<Forest> newForest(alternatePointer) ;
    
    newForest->InputForestElements(oriForest) ;
    newForest->SetBetweenMatListIndex(newBetweenMatListIndex) ;
    newForest->RebuildTrees(as<umat>(edgeMat) - 1) ;
    newForest->SetBetweenTransProbs(as<std::vector<mat>>(newBetweenTransProbs)) ; // We overwrite the between-cluster transition probabilities.
    newForest->InvalidateBetweenSolutions() ;
    newForest->ComputeLoglik() ;
    XPtr<Forest> p(newForest, false) ;
    
    return List::create(Named("logLik") = newForest->GetLoglik(), Named("solutionPointer") = p) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List newWithinTransProbsLogLik(SEXP ForestPointer, SEXP alternatePointer, List & newWithinTransProbs, IntegerMatrix & edgeMat, uint & numOpenMP, uint & newWithinMatListIndex) 
{
  //omp_set_num_threads(numOpenMP) ; 
  if (!(ForestPointer == NULL)) 
  {
    XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
    XPtr<Forest> newForest(alternatePointer) ;
    
    newForest->InputForestElements(oriForest) ;
    newForest->RebuildTrees(as<umat>(edgeMat) - 1) ;
    
    newForest->SetWithinTransProbs(as<std::vector<mat>>(newWithinTransProbs)) ;
    newForest->SetWithinMatListIndex(newWithinMatListIndex) ;
    newForest->InvalidateAllSolutions() ;
    newForest->ComputeLoglik() ;
    XPtr<Forest> p(newForest, false) ;
    return List::create(Named("logLik") = newForest->GetLoglik(), Named("solutionPointer") = p) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List withinClusNNIlogLik(SEXP ForestPointer, SEXP alternatePointer, IntegerMatrix & edgeMat, uint & MRCAofClusForNNI, uint & numMovesNNI, uint & numOpenMP) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointer == NULL)) 
  {
    XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
    XPtr<Forest> newForest(alternatePointer) ;
    
    newForest->InputForestElements(oriForest) ;
    newForest->RebuildTrees(as<umat>(edgeMat) - 1) ;
    
    AugTree * augTreePoint = newForest->GetForest().at(0) ;
    std::vector<uint> vertexIndexForNNI ;
    std::vector<uint> vertexIndexVec ;
    
    vertexIndexForNNI = augTreePoint->GetNNIverticesWithin(augTreePoint->GetVertexVector().at(MRCAofClusForNNI - 1)) ;
    
    for (uint counter = 0; counter < numMovesNNI; counter++)
    {
      unsigned long int rootForNNIindex = gsl_rng_uniform_int(newForest->GetRandomNumGenerator(), vertexIndexForNNI.size()) ;
      uvec placeholder(1) ;
      placeholder.at(0) = -1 ; // Probably a better way to do this than using a placeholder...
      vertexIndexVec = augTreePoint->GetTwoVerticesForNNI(newForest->GetRandomNumGenerator(), augTreePoint->GetVertexVector().at(vertexIndexForNNI.at(rootForNNIindex)), placeholder) ;
      newForest->RearrangeNNI(vertexIndexVec.at(0), vertexIndexVec.at(1)) ;
    }
    newForest->ComputeLoglik() ;
    umat newEdge = newForest->GetForest().at(0)->BuildEdgeMatrix(newForest->GetNumTips()) ; // All trees in the forest have the same hierarchy, hence the need to get the structure for only one of them.
    XPtr<Forest> p(newForest, false) ;
    return List::create(Named("logLik") = newForest->GetLoglik(),
                        Named("edge") = newEdge, Named("solutionPointer") = p) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List betweenClusNNIlogLik(SEXP ForestPointer, SEXP alternatePointer, NumericVector & clusterMRCAs, IntegerMatrix & edgeMat, uint & numMovesNNI, uint & numOpenMP) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointer == NULL)) 
  {
    XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
    XPtr<Forest> newForest(alternatePointer) ;
    
    newForest->InputForestElements(oriForest) ;
    newForest->RebuildTrees(as<umat>(edgeMat) - 1) ;
    
    AugTree * augTreePoint = newForest->GetForest().at(0) ;
    std::vector<uint> vertexIndexForNNI ;
    std::vector<uint> vertexIndexVec ;
    uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
    
    vertexIndexForNNI = augTreePoint->GetNNIverticesBetween(augTreePoint->GetVertexVector().at(newForest->GetNumTips()), clusterMRCAsRecast) ;
    
    for (uint counter = 0; counter < numMovesNNI; counter++)
    {
      unsigned long int rootForNNIindex = gsl_rng_uniform_int(newForest->GetRandomNumGenerator(), vertexIndexForNNI.size()) ;
      vertexIndexVec = augTreePoint->GetTwoVerticesForNNI(newForest->GetRandomNumGenerator(), augTreePoint->GetVertexVector().at(vertexIndexForNNI.at(rootForNNIindex)), clusterMRCAsRecast) ;
      
      for (auto & i : newForest->GetForest())
      {
        i->RearrangeTreeNNI(vertexIndexVec.at(0), vertexIndexVec.at(1), newForest->GetSolutionDictionary()) ;
      }
    }
    
    newForest->ComputeLoglik() ;
    umat newEdge = newForest->GetForest().at(0)->BuildEdgeMatrix(newForest->GetNumTips()) ;
    XPtr<Forest> p(newForest, false) ;
    return List::create(Named("logLik") = newForest->GetLoglik(),
                       Named("edge") = newEdge, Named("solutionPointer") = p) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List clusSplitMergeLogLik(SEXP ForestPointer, SEXP alternatePointer, IntegerVector & clusMRCAsToSplitOrMerge, List & withinTransProbsMats, List & betweenTransProbsMats, IntegerMatrix & edgeMat, uint & numOpenMP) 
{
  //omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointer == NULL))
  {
    XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
    uvec clusMRCAsToSplitOrMergeRecast = as<uvec>(clusMRCAsToSplitOrMerge) ;
    XPtr<Forest> newForest(alternatePointer) ;
    
    newForest->InputForestElements(oriForest) ;
    newForest->RebuildTrees(as<umat>(edgeMat) - 1) ;
    
    if (clusMRCAsToSplitOrMergeRecast.size() == 1) // Split move
    {
      std::vector<mat> betweenTransProbsMatsRecast = as<std::vector<mat>>(betweenTransProbsMats) ;
      newForest->HandleSplit(clusMRCAsToSplitOrMergeRecast.at(0)) ;
    }
    else
    {
      std::vector<mat> withinTransProbsMatsRecast = as<std::vector<mat>>(withinTransProbsMats) ;
      newForest->HandleMerge(clusMRCAsToSplitOrMergeRecast) ;
    }
    newForest->ComputeLoglik() ;
    XPtr<Forest> p(newForest, false) ;
    return List::create(Named("logLik") = newForest->GetLoglik(), Named("solutionPointer") = p) ;
  }
  else
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

void finalDeallocate(SEXP ForestPointer) // We need to explicitly deallocate the random number generator.
{
  XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
  gsl_rng_free(oriForest->GetRandomNumGenerator()) ;
  delete oriForest->GetSolutionDictionary() ;
  delete oriForest->GetAlignmentBinReference() ;
}

// [[Rcpp::export]]

void manualDeallocation(SEXP ForestPointer) // We need to explicitly deallocate the random number generator.
{
  XPtr<Forest> oriForest(ForestPointer) ; // Becomes a regular pointer again.
  Forest * pointerRecast = oriForest ;
  delete pointerRecast ;
}