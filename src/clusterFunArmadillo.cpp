// [[Rcpp::depends(RcppArmadillo)]]

//#include <gperftools/profiler.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <algorithm>
#include "AugTree.h"
#include <limits>


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

List logLikCpp(IntegerMatrix & edgeMat, NumericVector & clusterMRCAs, NumericVector & limProbsVec, List & withinTransMatList, List & betweenTransMatList, int numOpenMP, List alignmentBin, uint numTips, uint numLoci, SEXP pointerToForest)
{
  
  omp_set_num_threads(numOpenMP) ;
  if (!(pointerToForest == NULL)) // This syntax is really bad. Change once testing is complete!!!!
  {
    XPtr<Forest> Phylogenies(pointerToForest) ;
    Phylogenies->ComputeLoglik() ;
    return List::create(Named("logLik") = Phylogenies->GetLoglik(),
                        Named("ForestPointer") = Phylogenies) ;
  } 
  else
  {
    solutionDictionaryType solutionDictionary ;
    Forest Phylogenies(edgeMat, clusterMRCAs, alignmentBin, withinTransMatList, betweenTransMatList, limProbsVec, numTips, numLoci, solutionDictionary);
    Phylogenies.ComputeLoglik() ;
    XPtr<Forest> p(&Phylogenies, true) ;
    return List::create(Named("logLik") = Phylogenies.GetLoglik(),
                        Named("ForestPointer") = p) ;
  }
  
  // std::vector<std::vector<vec>>* mySolutions = new std::vector<std::vector<vec>> ;
  // for (auto & aTree : Phylogenies.GetForest()) 
  // {
  //   mySolutions->push_back(aTree->GetSolutionsFromTree()) ;
  // }
  // 
  // Rcpp::XPtr<std::vector<std::vector<vec>>> p(mySolutions, true) ;
  //Rcpp::XPtr<Forest> p(&Phylogenies, true) ; 
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

SEXP newBetweenTransProbsLogLik(SEXP ForestPointer, List & newBetweenTransProbs) 
{
  if (!(ForestPointer == NULL)) 
  {
    XPtr<Forest> Phylogenies(ForestPointer) ; // Becomes a regular pointer again.
    std::vector<mat> newBetweenTransProbsRecast = as<std::vector<mat>>(newBetweenTransProbs) ;
    Phylogenies->AmendBetweenTransProbs(newBetweenTransProbsRecast) ;
    Phylogenies->ComputeLoglik() ;
    return List::create(Named("logLik") = Phylogenies->GetLoglik(),
                        Named("ForestPointer") = ForestPointer) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

SEXP newWithinTransProbsLogLik(SEXP ForestPointer, List newWithinTransProbs, IntegerVector clusterMRCAs) 
{
  if (!(ForestPointer == NULL)) 
  {
    XPtr<Forest> Phylogenies(ForestPointer) ; // Becomes a regular pointer again.
    std::vector<mat> newWithinTransProbsRecast = as<std::vector<mat>>(newWithinTransProbs) ;
    uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
    Phylogenies->AmendWithinTransProbs(newWithinTransProbsRecast, clusterMRCAsRecast) ;
    Phylogenies->ComputeLoglik() ;
    return List::create(Named("logLik") = Phylogenies->GetLoglik(),
                        Named("ForestPointer") = ForestPointer) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

SEXP withinClusNNIlogLik(SEXP ForestPointer, uint MRCAofClusForNNI) 
{
  // TO_DO
}

// [[Rcpp::export]]

SEXP betweenClusNNIlogLik(SEXP ForestPointer) 
{
  // TO_DO
}

// [[Rcpp::export]]

SEXP clusSplitMergeLogLik(SEXP ForestPointer, IntegerVector clusMRCAsToSplitOrMerge, List withinTransProbsMats, List betweenTransProbsMats) 
{
  XPtr<Forest> Phylogenies(ForestPointer) ; // Becomes a regular pointer again.
  if (!(ForestPointer == NULL)) 
  {
    uvec clusMRCAsToSplitOrMergeRecast = as<uvec>(clusMRCAsToSplitOrMerge) ;
    
    if (clusMRCAsToSplitOrMergeRecast.size() == 1) // Split move
    {
      std::vector<mat> betweenTransProbsMatsRecast = as<std::vector<mat>>(betweenTransProbsMats) ;
      Phylogenies->HandleSplit(clusMRCAsToSplitOrMergeRecast.at(0), betweenTransProbsMatsRecast) ;
    }
    else
    {
      std::vector<mat> withinTransProbsMatsRecast = as<std::vector<mat>>(withinTransProbsMats) ;
      Phylogenies->HandleMerge(clusMRCAsToSplitOrMergeRecast, withinTransProbsMatsRecast) ;
    }
    Phylogenies->ComputeLoglik() ;
    return List::create(Named("logLik") = Phylogenies->GetLoglik(),
                        Named("ForestPointer") = ForestPointer) ;
  }
  else
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
  
}