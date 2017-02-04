#include <RcppArmadillo.h>
#include "PhyloAlpha.h"
using namespace Rcpp;
using namespace arma;

phylogenyAlpha::phylogenyAlpha(const NumericMatrix & edgeMat, const CharacterMatrix & alignmentAlphaMat, const NumericVector & logLimProbsVec, const List & logTransMatList, const int numOpenMP, const CharacterVector & equivVec, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) : phylo(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsVec) {
  
  alignmentRcpp = alignmentAlphaMat ;
  equivalency = Rcpp::as<std::vector<std::string>>(equivVec) ;
  convertSTL() ;
  defineMap() ;
  convertNucleoToNum() ;
}

// This is the faster constructor, that does not convert the character alignment into matrices of identity vectors, because it's provided by the user (alignmentList).
phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::List & alignmentList, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numOpenMP, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) : phylo(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsVec) {
  
  alignmentBin = alignmentList ;
  convertAlignmentList(alignmentList) ; 
  
  if (!internalFlag) 
  {
    numTips = as<NumericMatrix>(alignmentBin[0]).ncol() ;
    numUniqueLoci = alignmentBin.size() ;
  } 
  else
  { // Tip vectors depend on the rate category.
    numTips = alignmentBin.size() ;
    numUniqueLoci = as<NumericMatrix>(as<List>(alignmentBin[0])[0]).ncol() ;
  }
  Nnode = edgeMat.nrow() - numTips + 1 ;
  compUpdateVec() ;
}

phylogenyAlpha::phylogenyAlpha():phylo() {}

void initializeGraph() { // Note that this function assumes that the ordering of elements in alignmentBin matches the tip numbers in edge!
  //TO_DO
}

void phylogenyAlpha::convertSTL() {
  
  alignmentAlpha.resize(alignmentRcpp.ncol()) ;
  
  for (int i = 0; i < alignmentRcpp.ncol(); i++) {
    
    Rcpp::CharacterVector locusCol = alignmentRcpp(Rcpp::_, i) ;
    alignmentAlpha[i] = Rcpp::as<std::vector<std::string>>(locusCol) ;
    
  }
}

void phylogenyAlpha::defineMap() {
  
  for (uint i = 0; i < numStates; i++) {
    
    vec unitVec(numStates, fill::zeros) ;
    unitVec(i) = 1 ;
    std::string label = equivalency[i] ;
    mymap[label] = unitVec ;
  }
  if (numStates == 4) 
  { // If we have four alphanumeric states, the assumption is that we're looking at a nucleotide alignment with states "a", "t", "c", "g".
    // Now, we handle ambiguities.
    mymap["r"] = (mymap["a"] + mymap["g"]) ;
    mymap["y"] = (mymap["c"] + mymap["t"]) ;
    mymap["s"] = (mymap["c"] + mymap["g"]) ;
    mymap["w"] = (mymap["a"] + mymap["t"]) ;
    mymap["k"] = (mymap["g"] + mymap["t"]) ;
    mymap["m"] = (mymap["a"] + mymap["c"]) ;
    mymap["b"] = (mymap["c"] + mymap["g"] + mymap["t"]) ;
    mymap["d"] = (mymap["a"] + mymap["g"] + mymap["t"]) ;
    mymap["h"] = (mymap["a"] + mymap["c"] + mymap["t"]) ;
    mymap["v"] = (mymap["a"] + mymap["c"] + mymap["g"]) ;
    mymap["-"] = (mymap["a"] + mymap["t"] + mymap["c"] + mymap["g"]) ;
    mymap["n"] = (mymap["a"] + mymap["t"] + mymap["c"] + mymap["g"]) ;
    mymap["."] = (mymap["a"] + mymap["t"] + mymap["c"] + mymap["g"]) ;
  } 
}

void phylogenyAlpha::convertNucleoToNum()  {
  
  List alignmentBinList(numLoci) ;
  std::vector<std::string> locusCol ;
  mat containerMat(numStates, alignmentRcpp.nrow(), fill::zeros) ;
  for (uint locusIndex = 0 ; locusIndex < numLoci ; locusIndex++) {
    
    locusCol = alignmentAlpha[locusIndex] ;
    for (uint i = 0; i < containerMat.n_cols; i++) {
      
      containerMat.col(i) = mymap[locusCol[i]] ;
    }
    alignmentBinList(locusIndex) = as<NumericMatrix>(wrap(containerMat)) ;
  }
  
  alignmentBin = alignmentBinList ; 
}
