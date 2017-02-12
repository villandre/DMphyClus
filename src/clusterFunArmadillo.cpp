// [[Rcpp::depends(RcppArmadillo)]]

//#include <gperftools/profiler.h>
#include <RcppArmadillo.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <algorithm>
#include "AugTree.h"

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

List logLikCpp(IntegerMatrix & edgeMat, NumericVector & clusterMRCAs, NumericVector & limProbsVec, List & withinTransMatList, List & betweenTransMatList, int numOpenMP, List alignmentBin, int numTips)
{
  omp_set_num_threads(numOpenMP) ;
  solutionDictionaryType solutionDictionary ;
  Forest Phylogenies(edgeMat, clusterMRCAs, alignmentBin, withinTransMatList, betweenTransMatList, limProbsVec, numTips, solutionDictionary);
  Phylogenies.ComputeLoglik() ; 
  return List::create(Named("logLik") = Phylogenies.GetLoglik(),
                      Named("intermediateNodes") = wrap(0)) ;
}

std::unordered_map<std::string, Col<uvec>> defineMap(std::vector<std::string> & equivalency) 
{
  uint numStates = equivalency.size() ;
  std::unordered_map<std::string, Col<uvec>> myMap ;
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

SEXP getConvertedAlignment(int numOpenMP, SEXP & equivVector, CharacterMatrix & alignmentAlphaMat)
{
  std::vector<std::vector<uvec>> alignmentBin ;
  alignmentBin.reserve(alignmentAlphaMat.ncol()) ;
  std::vector<std::string> equivVectorAlpha = as<std::vector<std::string>>(equivVector) ;
  std::unordered_map<std::string, Col<uvec>> siteMap = defineMap(equivVectorAlpha) ;
  
  for (uint i = 0; i < alignmentBin.size(); i++) 
  {
    alignmentBin.at(i).reserve(alignmentAlphaMat.nrow()) ;
    std::transform(alignmentAlphaMat.begin()+i*alignmentAlphaMat.nrow(),
                   alignmentAlphaMat.begin()+(i+1)*alignmentAlphaMat.nrow(),
                   alignmentBin.at(i).begin(),
                   [& siteMap] (std::string & myLetter) { return siteMap[myLetter] ;}) ;
  }
  return wrap(alignmentBin) ;
}
