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

List logLikCpp(IntegerMatrix & edgeMat, NumericVector & clusterMRCAs, NumericVector & limProbsVec, List & withinTransMatList, List & betweenTransMatList, int numOpenMP, List alignmentBin, uint numTips, uint numLoci)
{
  omp_set_num_threads(numOpenMP) ;
  
  solutionDictionaryType solutionDictionary ;
  std::vector<Forest *>* PhyloForestsPoint = new std::vector<Forest *>;
  PhyloForestsPoint->reserve(2) ;
  Forest * PhylogeniesPoint1 = new Forest(edgeMat, clusterMRCAs, alignmentBin, withinTransMatList, betweenTransMatList, limProbsVec, numTips, numLoci, solutionDictionary);
  PhylogeniesPoint1->ComputeLoglik() ;
  PhyloForestsPoint->push_back(PhylogeniesPoint1) ;
  Forest * PhylogeniesPoint2 = new Forest(edgeMat, clusterMRCAs, alignmentBin, withinTransMatList, betweenTransMatList, limProbsVec, numTips, numLoci, solutionDictionary);
  *PhylogeniesPoint2 = *PhylogeniesPoint1 ;
  PhyloForestsPoint->push_back(PhylogeniesPoint2) ;
  cout << "First pointer: " << PhyloForestsPoint->at(0) << " Second pointer: " << PhyloForestsPoint->at(1) << "\n" ;
  XPtr<std::vector<Forest *>> p(PhyloForestsPoint, true) ;
  cout << "Done! \n" ;
  return List::create(Named("logLik") = PhylogeniesPoint1->GetLoglik(),
                      Named("ForestPointer") = p) ;
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

List newBetweenTransProbsLogLik(SEXP ForestPointerVec, List & newBetweenTransProbs, int numOpenMP) 
{
  omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointerVec == NULL)) 
  {
    XPtr<std::vector<Forest *>> TwoForests(ForestPointerVec) ; // Becomes a regular pointer again.
    cout << "First pointer: " << TwoForests->at(0) << " Second pointer: " << TwoForests->at(1) << "\n" ;
    TwoForests->at(0)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("Before:") ;
    std::vector<mat> newBetweenTransProbsRecast = as<std::vector<mat>>(newBetweenTransProbs) ;
    TwoForests->at(1)->AmendBetweenTransProbs(newBetweenTransProbsRecast) ;
    TwoForests->at(1)->ComputeLoglik() ;
    TwoForests->at(0)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("After:") ;
    TwoForests->at(1)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("After at new loc.:") ;
    
    return List::create(Named("logLik") = TwoForests->at(1)->GetLoglik()) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List newWithinTransProbsLogLik(SEXP ForestPointerVec, List newWithinTransProbs, IntegerVector clusterMRCAs, int numOpenMP) 
{
  if (!(ForestPointerVec == NULL)) 
  {
    XPtr<std::vector<Forest *>> TwoForests(ForestPointerVec) ; // Becomes a regular pointer again.
    std::vector<mat> newWithinTransProbsRecast = as<std::vector<mat>>(newWithinTransProbs) ;
    uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
    TwoForests->at(1)->AmendWithinTransProbs(newWithinTransProbsRecast, clusterMRCAsRecast) ;
    TwoForests->at(1)->ComputeLoglik() ;
    return List::create(Named("logLik") = TwoForests->at(1)->GetLoglik()) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List withinClusNNIlogLik(SEXP ForestPointerVec, uint MRCAofClusForNNI, uint numMovesNNI, int numOpenMP) 
{
  omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointerVec == NULL)) 
  {
    XPtr<std::vector<Forest *>> TwoForests(ForestPointerVec) ; // Becomes a regular pointer again.
    AugTree * augTreePoint = TwoForests->at(1)->GetForest().at(0) ;
    std::vector<uint> vertexIndexForNNI ;
    std::vector<uint> vertexIndexVec ;
    
    vertexIndexForNNI = TwoForests->at(1)->GetForest().at(0)->GetNNIverticesWithin(augTreePoint->GetVertexVector().at(MRCAofClusForNNI - 1)) ;
    
    for (uint counter = 0; counter < numMovesNNI; counter++)
    {
      unsigned long int rootForNNIindex = gsl_rng_uniform_int(TwoForests->at(1)->GetRandomNumGenerator(), vertexIndexForNNI.size()) ;
      uvec placeholder(1) ;
      placeholder.at(0) = -1 ; // Probably a better way to do this than using a placeholder...
      vertexIndexVec = augTreePoint->GetTwoVerticesForNNI(TwoForests->at(1)->GetRandomNumGenerator(), augTreePoint->GetVertexVector().at(vertexIndexForNNI.at(rootForNNIindex)), placeholder) ;
      
      for (auto & i : TwoForests->at(1)->GetForest())
      {
        i->RearrangeTreeNNI(vertexIndexVec.at(0), vertexIndexVec.at(1)) ;
      }
    }
    TwoForests->at(1)->ComputeLoglik() ;
    umat newEdge = TwoForests->at(1)->GetForest().at(0)->BuildEdgeMatrix() ; // All trees in the forest have the same hierarchy, hence the need to get the structure for only one of them.
    return List::create(Named("logLik") = TwoForests->at(1)->GetLoglik(),
                        Named("edge") = newEdge) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List betweenClusNNIlogLik(SEXP ForestPointerVec, uint numMovesNNI, int numOpenMP, NumericVector & clusterMRCAs) 
{
  omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointerVec == NULL)) 
  {
    XPtr<std::vector<Forest *>> TwoForests(ForestPointerVec) ; // Becomes a regular pointer again.
    AugTree * augTreePoint = TwoForests->at(1)->GetForest().at(0) ;
    std::vector<uint> vertexIndexForNNI ;
    std::vector<uint> vertexIndexVec ;
    uint numTips = augTreePoint->GetNumTips() ;
    uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
    
    vertexIndexForNNI = TwoForests->at(1)->GetForest().at(0)->GetNNIverticesBetween(augTreePoint->GetVertexVector().at(numTips), clusterMRCAsRecast) ;
    
    for (uint counter = 0; counter < numMovesNNI; counter++)
    {
      unsigned long int rootForNNIindex = gsl_rng_uniform_int(TwoForests->at(1)->GetRandomNumGenerator(), vertexIndexForNNI.size()) ;
      vertexIndexVec = augTreePoint->GetTwoVerticesForNNI(TwoForests->at(1)->GetRandomNumGenerator(), augTreePoint->GetVertexVector().at(vertexIndexForNNI.at(rootForNNIindex)), clusterMRCAsRecast) ;
      
      for (auto & i : TwoForests->at(1)->GetForest())
      {
        i->RearrangeTreeNNI(vertexIndexVec.at(0), vertexIndexVec.at(1)) ;
      }
    }
    TwoForests->at(1)->ComputeLoglik() ;
    umat newEdge = TwoForests->at(1)->GetForest().at(0)->BuildEdgeMatrix() ;
    return List::create(Named("logLik") = TwoForests->at(1)->GetLoglik(),
                       Named("edge") = newEdge) ;
  } 
  else 
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

List clusSplitMergeLogLik(SEXP ForestPointerVec, IntegerVector & clusMRCAsToSplitOrMerge, List & withinTransProbsMats, List & betweenTransProbsMats, int numOpenMP) 
{
  omp_set_num_threads(numOpenMP) ;
  if (!(ForestPointerVec == NULL))
  {
    XPtr<std::vector<Forest *>> TwoForests(ForestPointerVec) ; // Becomes a regular pointer again.
    uvec clusMRCAsToSplitOrMergeRecast = as<uvec>(clusMRCAsToSplitOrMerge) ;
    
    if (clusMRCAsToSplitOrMergeRecast.size() == 1) // Split move
    {
      std::vector<mat> betweenTransProbsMatsRecast = as<std::vector<mat>>(betweenTransProbsMats) ;
      TwoForests->at(1)->HandleSplit(clusMRCAsToSplitOrMergeRecast.at(0), betweenTransProbsMatsRecast) ;
    }
    else
    {
      std::vector<mat> withinTransProbsMatsRecast = as<std::vector<mat>>(withinTransProbsMats) ;
      TwoForests->at(1)->HandleMerge(clusMRCAsToSplitOrMergeRecast, withinTransProbsMatsRecast) ;
    }
    TwoForests->at(1)->ComputeLoglik() ;
    return List::create(Named("logLik") = TwoForests->at(1)->GetLoglik()) ;
  }
  else
  {
    throw ::Rcpp::exception("pointer is null." ) ;
  }
}

// [[Rcpp::export]]

void copyForestElements(bool keepOld, SEXP ForestVecPointer)
{
  XPtr<std::vector<Forest *>> PhyloForestsPoint(ForestVecPointer) ;
  if (keepOld)
  {
    cout << "Keeping old... \n" ;
    cout << "Before: \n" ;
    PhyloForestsPoint->at(0)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("TransMatrix:") ;
    PhyloForestsPoint->at(1)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("TransMatrix:") ;
    //std::copy(PhyloForestsPoint->begin(), PhyloForestsPoint->begin(), PhyloForestsPoint->begin()+1) ;
    *PhyloForestsPoint->at(1)=*PhyloForestsPoint->at(0) ;
    cout << "After: \n" ;
    PhyloForestsPoint->at(0)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("TransMatrix:") ;
    PhyloForestsPoint->at(1)->GetForest().at(0)->GetVertexVector().at(6)->GetTransMatrix().print("TransMatrix:") ;
  }
  else
  {
    *PhyloForestsPoint->at(0)=*PhyloForestsPoint->at(1) ;
  }
}