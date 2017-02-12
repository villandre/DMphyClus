// [[Rcpp::depends(RcppArmadillo)]]

//#include <gperftools/profiler.h>
#include <RcppArmadillo.h>
#include <vector>
#include <stdexcept>
#include <omp.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <sys/resource.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#include "PhyloAlpha.h"
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

// // [[Rcpp::export]]a
// 
// SEXP logLikCppToWrap(NumericMatrix & edgeMat, NumericVector & logLimProbsVec, List & logTransMatList, int numOpenMP, SEXP & equivVector, List alignmentBin, const bool returnMatIndic, const bool internalFlag, const NumericVector sitePatterns) {
// 
//   omp_set_num_threads(numOpenMP) ;
// 
//   uint numLoci ;
//   uvec sitePatternsAmended ;
//   if (internalFlag) 
//   { 
//     numLoci = as<NumericMatrix>(as<List>(alignmentBin[0])[0]).ncol() ;
//     sitePatternsAmended.resize(numLoci) ;
//     sitePatternsAmended = conv_to<uvec>::from(linspace(1, numLoci, numLoci))  ;
//   } 
//   else
//   {
//     sitePatternsAmended = as<uvec>(sitePatterns) ;
//   }
//   SEXP container ;
//   
//   //ProfilerStart("/home/villandre/profileOut.out") ;
//   phylogenyAlpha phyloObject(edgeMat, alignmentBin, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsAmended) ;
//   phyloObject.logLikPhylo(returnMatIndic);
//   //ProfilerStop() ;
//   if (returnMatIndic) 
//   {
//     container = phyloObject.getPruningMat();
//   } 
//   else 
//   {
//     container = phyloObject.getLogLik();
//   }
//   return container ;
// }

// [[Rcpp::export]]

SEXP getConvertedAlignmentToWrap(int numOpenMP, SEXP & equivVector, Rcpp::CharacterMatrix & alignmentAlphaMat, NumericVector & sitePatterns) {
  
  double src[] = {3, 3, 1, 2};
  NumericMatrix placeholderMat(2, 2, src) ; //src is considered an iterator.
  NumericVector placeholderVec(as<CharacterVector>(equivVector).size(),1000.0) ;
  List placeholderList(1) ; // This is a placeholder: 1 is an arbitrary number.
  placeholderList[0] = Rcpp::NumericMatrix(2,2) ;
  bool returnMatIndic = false ;
  bool internalFlag = false ;
  
  phylogenyAlpha phyloObject(placeholderMat, alignmentAlphaMat, placeholderVec, placeholderList, numOpenMP, equivVector, returnMatIndic, internalFlag, as<uvec>(sitePatterns)) ;
  return phyloObject.getAlignmentBin() ;
}

template<class... Conts>
auto zip_range(Conts&... conts)
  -> decltype(boost::make_iterator_range(
      boost::make_zip_iterator(boost::make_tuple(conts.begin()...)),
      boost::make_zip_iterator(boost::make_tuple(conts.end()...))))
  {
    return {boost::make_zip_iterator(boost::make_tuple(conts.begin()...)),
            boost::make_zip_iterator(boost::make_tuple(conts.end()...))};
  }

// [[Rcpp::export]]

List getSitePatterns(List alignmentBin) {
  
  std::size_t myHash ;
  std::vector<std::size_t> hashVector ;
  NumericVector sitePatterns(alignmentBin.size()) ;
  std::vector<NumericMatrix> uniqueDNAdataBin ;
  auto hashIterator = hashVector.begin() ;
  std::vector<NumericMatrix> convertedAlignmentBin = as<std::vector<NumericMatrix> >(alignmentBin) ;
  
  for(auto&& i : zip_range(convertedAlignmentBin, sitePatterns)) {
    
    myHash = boost::hash_range(as<NumericMatrix>(i.get<0>()).begin(), as<NumericMatrix>(i.get<0>()).end()) ;
    
    hashIterator = std::find(hashVector.begin(), hashVector.end(), myHash) ;
    i.get<1>() = std::distance(hashVector.begin(), hashIterator) + 1 ;
    
    if (hashIterator == hashVector.end()) 
    {
      hashVector.push_back(myHash) ;
      uniqueDNAdataBin.push_back(as<NumericMatrix>(i.get<0>())) ;
    }
  }
  
  return List::create(Named("uniqueDNAdataBin") = wrap(uniqueDNAdataBin),
                      Named("sitePatterns") = wrap(sitePatterns));
}

// [[Rcpp::export]]

List logLikCppToWrap(IntegerMatrix & edgeMat, NumericVector & clusterMRCAs, NumericVector & limProbsVec, List & withinTransMatList, List & betweenTransMatList, int numOpenMP, List alignmentBin, int numTips)
{
  omp_set_num_threads(numOpenMP) ;
  solutionDictionaryType solutionDictionary ;
  Forest Phylogenies(edgeMat, clusterMRCAs, alignmentBin, withinTransMatList, betweenTransMatList, limProbsVec, numTips, solutionDictionary);
  Phylogenies.ComputeLoglik() ; 
  return List::create(Named("logLik") = Phylogenies.GetLoglik(),
                      Named("intermediateNodes") = wrap(Phylogenies.GetSolutions())) ;
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
