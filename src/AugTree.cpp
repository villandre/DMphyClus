#include "AugTree.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include<vector.h>
#include<memory.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#include <unordered_map>

using namespace Rcpp ;
using namespace arma ;

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}


AugTree::AugTree(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, const mat & withinTransProbMatrix, const mat & betweenTransProbMatrix, const std::vector<uvec> & alignmentBinOneLocus, const NumericVector & limProbs, const uint numTips, const uint rateCategIndex) 
{
  _numTips = numTips ;
  _limProbs = as<Col<long double>>(limProbs) ;
  umat edgeMatrixCopy(as<umat>(edgeMatrix)) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
  InitializeTips(alignmentBinOneLocus) ;
  AssociateTransProbMatrices(clusterMRCAs, withinTransProbMatrix, betweenTransProbMatrix, rateCategIndex) ;
  IdentifyPatterns(_tree.at(numTips)) ;
}

void AugTree::AssociateTransProbMatrices(const NumericVector & clusterMRCAs, const mat & withinTransProbMatrix, const mat & betweenTransProbMatrix, const uint rateCategIndex) {
  
  // By default all nodes are considered between clusters.
  for (auto & i : _tree) 
  {
    i->SetTransProbMatrix(betweenTransProbMatrix, rateCategIndex, false) ;
  }
  for (auto & i : clusterMRCAs) 
  { 
    if (i > _numTips) { // Again, clusterMRCAs is based on the R convention, hence >, and not >=.
      for (auto & j : _tree.at(i-1)->GetChildren()) // The -1 accounts for the difference between R and C++ in terms of indices.
      {
        BindMatrixChildren(j, withinTransProbMatrix, rateCategIndex, true) ; 
      }
    }
  }
}

void AugTree::BindMatrixChildren(TreeNode * vertex, const mat & transProbMatrix, const uint rateCategIndex, const bool withinCluster)
{
  vertex->SetTransProbMatrix(transProbMatrix, rateCategIndex, withinCluster) ;
  if (!(vertex->GetChildren()[0] == NULL)) { // A null pointer indicates that we've reached an input node.
    for (auto & i : vertex->GetChildren()) 
    {
      BindMatrixChildren(i, transProbMatrix, rateCategIndex, withinCluster) ; 
    }
  } 
}
void AugTree::BuildTree(umat & edgeMatrix) 
{
  _tree.reserve(edgeMatrix.n_rows + 1) ;
  
  // We create the tips. Note that tip 1 should correspond to vertex 1 in the original (the one in the phylo object) edgeMatrix
  
  for (uint i = 0; i < _numTips; i++) {
    InputNode newNode{} ;  
    _tree.push_back(&newNode) ;
  } ;
  
  // We add the internal nodes
  for (uint i = 0 ; i < edgeMatrix.n_rows - _numTips + 1; i++) {
    IntermediateNode myNode{} ;
    _tree.at(i) = &myNode ;
  } ;
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _tree.size(); i++) {
    _tree.at(i)->SetId(i+1) ; // R want edge labels that start at 1.
  } ;
  
  // The vertices are all disjoint, the next loop defines their relationships
  for (umat::iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end(); iter = iter + 2) 
  { 
    _tree.at(*iter)->AddChild(_tree.at(*(iter+1))) ;
    _tree.at(*(iter+1))->SetParent(_tree.at(*iter)) ;
  }
}

SEXP AugTree::BuildEdgeMatrix() 
{
  //TO_DO
}

void AugTree::IdentifyPatterns(TreeNode * vertex) 
{
  if (vertex->GetChildren()[0] != NULL) { 
    for (auto & i : vertex->GetChildren())
    {
      if (i->CanSolve())
      {
        i->SetPattern() ;
        i->ToggleSolved() ;
      }
      else
      {
        IdentifyPatterns(i) ;
      }
    }
  }
  // The nodes had been marked as solved only so that we could get a pattern indicator for each. They are not really solved and so, must be marked as such.
  // Note that, as it should, ToggleSolved() does nothing for InputNodes.
  for (auto & i : _tree)
  {
    i->ToggleSolved() ;
  }
}

void AugTree::SolveRoot(solutionDictionaryType & solutionDictionary) {
  PatternLookup(solutionDictionary, _tree.at(_numTips - 1)) ;
  if (!(_tree.at(_numTips - 1)->IsSolved()))
  {
    if (!_tree.at(_numTips - 1)->CanSolve())
    {
      for (auto & i : _tree.at(_numTips - 1)->GetChildren())
      {
        TrySolve(i) ;
      }
    }
    else
    {
      _tree.at(_numTips - 1)->ComputeSolution() ;
    }
  }
  _likelihood = dot(_tree.at(_numTips - 1)->GetSolution(), _limProbs) ; 
}

void AugTree::SetPatterns() {
  //TO_DO
}

void AugTree::InitializeFromDictionary() 
{
  //TO_DO
}

void AugTree::InitializeTips(const std::vector<uvec> & alignmentBinOneLocus) 
{
  std::vector<TreeNode *>::iterator treeIter = _tree.begin() ;
  for (auto & i : alignmentBinOneLocus) 
  {
    (*treeIter)->SetInput(i) ; 
    treeIter++ ;
  }
}

void AugTree::PatternLookup(solutionDictionaryType & solutionDictionary, TreeNode * currentNode) {
  if (!currentNode->IsSolved()) { // If the node is already solved, no need to update it with a stored pattern.
    solutionDictionaryType::iterator patternIter=solutionDictionary.find(currentNode->GetDictionaryKey()) ;
    if (patternIter == solutionDictionary.end()) 
    {
      for (auto & i : currentNode->GetChildren())
      {
        PatternLookup(solutionDictionary, i) ;
      }
    } 
    else 
    {
      currentNode->SetSolution((*patternIter).second) ;
      currentNode->ToggleSolved() ;
    }
  }
}

void TrySolve(TreeNode * vertex)  
{
  if (!(vertex->IsSolved())) 
  {
    if (vertex->CanSolve())
    {
      vertex->ComputeSolution() ;
    } 
    else
    {
      for (auto & i : vertex->GetChildren())
      {
        TrySolve(i) ;
      }
    }
  }
}

Forest::Forest(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, const List & alignmentBin, const List & withinTransProbMatList, const List & betweenTransProbMatList, const NumericVector & limProbs, const uint numTips) 
{
  _numLoci = alignmentBin.size() ;
  _numRateCats = withinTransProbMatList.size() ; 
  std::vector<AugTree> phyloForest ;
  phyloForest.reserve(alignmentBin.size()*withinTransProbMatList.size()) ; // We will have one tree per locus per rate category.
  
  for (auto & i : alignmentBin) 
  {
    for (auto j : zip(withinTransProbMatList, betweenTransProbMatList)) 
    {
      phyloForest.push_back(AugTree(edgeMatrix, clusterMRCAs, j.get<0>(), j.get<1>(), as<std::vector<uvec>>(i), limProbs, numTips)) ;
    }
  } // In the forest, elements 0,..., numRateCats - 1 are for locus 1, elements numRateCats,..., 2*numRateCats - 1 are for locus 2, and so on.
  _forest = phyloForest ;
}

void Forest::ComputeLoglik() 
{
  for (auto & forestIter : _forest) 
  {
    forestIter.SolveRoot(_solutionDictionary) ;
  }
  // Now, we must average likelihoods across rate categories for each locus, log the output, and sum the resulting logs.
  Col<long double> rateAveragedLogLiks(_numLoci) ;
  Col<long double> likAcrossRatesLoci(_forest.size()) ;
  std::transform(_forest.begin(), _forest.end(), likAcrossRatesLoci.begin(), [] (const AugTree & myTree) {return myTree.GetLikelihood() ;}) ;
  for (uint i = 0; i < rateAveragedLogLiks.size(); i++) {
    rateAveragedLogLiks.at(i) = log(mean(likAcrossRatesLoci.rows(i, i + _numRateCats - 1))) ;
  }
  _loglik = sum(rateAveragedLogLiks) ;
}

void Forest::NNmovePropagate() 
{
  
}



