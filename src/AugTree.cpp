#include "AugTree.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

using namespace Rcpp ;
using namespace arma ;

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}


AugTree::AugTree(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, const mat & withinTransProbMatrix, const mat & betweenTransProbMatrix, const std::vector<uvec> & alignmentBinOneLocus, const NumericVector & limProbs, const uint numTips, const uint rateCategIndex, solutionDictionaryType & solutionDictionary) 
{
  _numTips = numTips ;
  _limProbs = as<Col<long double>>(limProbs) ;
  _rateCateg = rateCategIndex ;
  umat edgeMatrixCopy(as<umat>(edgeMatrix)) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
  InitializeTips(alignmentBinOneLocus) ;
  AssociateTransProbMatrices(clusterMRCAs, withinTransProbMatrix, betweenTransProbMatrix) ;
  ComputeKeys(_tree.at(numTips-1), solutionDictionary) ; // We start obtaining keys at the root.
}

void AugTree::AssociateTransProbMatrices(const NumericVector & clusterMRCAs, const mat & withinTransProbMatrix, const mat & betweenTransProbMatrix) {
  
  // By default all nodes are considered between clusters.
  for (auto & i : _tree) 
  {
    i->SetTransProbMatrix(betweenTransProbMatrix, _rateCateg, false) ;
  }
  for (auto & i : clusterMRCAs) 
  { 
    if (i > _numTips) { // Again, clusterMRCAs is based on the R convention, hence >, and not >=.
      for (auto & j : _tree.at(i-1)->GetChildren()) // The -1 accounts for the difference between R and C++ in terms of indices.
      {
        BindMatrixChildren(j, withinTransProbMatrix, true) ; 
      }
    }
  }
}

void AugTree::BindMatrixChildren(TreeNode * vertex, const mat & transProbMatrix, const bool withinCluster) 
{
  vertex->SetTransProbMatrix(transProbMatrix, _rateCateg, withinCluster) ;
  if (!(vertex->GetChildren()[0] == NULL)) { // A null pointer indicates that we've reached an input node.
    for (auto & i : vertex->GetChildren()) 
    {
      BindMatrixChildren(i, transProbMatrix, withinCluster) ; 
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
  return wrap(0) ;
}

void AugTree::ComputeKeys(TreeNode * vertex, solutionDictionaryType & solutionDictionary) 
{
  if (vertex->CanSolve()) // If a vertex can be solved, then its children have patterns assigned to them.
  {
    vertex->DeriveKey(solutionDictionary) ;
    vertex->ToggleSolved() ;
  } 
  else // Vertex cannot be a tip, else, CanSolve would have returned true.
  {
    for (auto & i : vertex->GetChildren())
    {
      ComputeKeys(i, solutionDictionary) ;
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
  if (!currentNode->IsSolved()) 
    { // If the node is already solved, no need to update it with a stored pattern.
    solutionDictionaryType::iterator patternIter = solutionDictionary.find(currentNode->GetDictionaryKey()) ;
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
  if (!(vertex->IsSolved())) // Could be solved because of the pattern lookup.
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

Forest::Forest(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, const List & alignmentBin, const List & withinTransProbMatList, const List & betweenTransProbMatList, const NumericVector & limProbs, const uint numTips, solutionDictionaryType & solutionDictionary)
{
  _numLoci = alignmentBin.size() ;
  _numRateCats = withinTransProbMatList.size() ;
  _solutionDictionary = solutionDictionary ;
  _forest.reserve(alignmentBin.size()*withinTransProbMatList.size()) ;
  
  for (auto & i : alignmentBin) // Iterating on loci...
  { 
    uint rateCategIndex = 0 ;
    for (auto j : zip(withinTransProbMatList, betweenTransProbMatList)) 
    {
      _forest.push_back(AugTree(edgeMatrix, clusterMRCAs, j.get<0>(), j.get<1>(), as<std::vector<uvec>>(i), limProbs, numTips, rateCategIndex, _solutionDictionary)) ;
      rateCategIndex++ ;
    }
  } // In the forest, elements 0,..., numRateCats - 1 are for locus 1, elements numRateCats,..., 2*numRateCats - 1 are for locus 2, and so on.
}

void Forest::ComputeLoglik() 
{
  #pragma omp parallel for
  for (std::vector<AugTree>::iterator forestIter = _forest.begin(); forestIter < _forest.end(); forestIter++) // This syntax is compatible with openMP, unlike the more conventional 'for (auto & i : myVec')
  {
    forestIter->SolveRoot(_solutionDictionary) ;
  }
  // Now, we must average likelihoods across rate categories for each locus, log the output, and sum the resulting logs.
  Col<long double> rateAveragedLogLiks(_numLoci) ;
  Col<long double> likAcrossRatesLoci(_forest.size()) ;
  std::transform(_forest.begin(), _forest.end(), likAcrossRatesLoci.begin(), [] (const AugTree & myTree) {return myTree.GetLikelihood() ;}) ;
  for (uint i = 0; i < rateAveragedLogLiks.size(); i++) 
  {
    rateAveragedLogLiks.at(i) = log(mean(likAcrossRatesLoci.rows(i, i + _numRateCats - 1))) ;
  }
  _loglik = sum(rateAveragedLogLiks) ;
}

void Forest::NNmovePropagate() 
{
 //TO_DO 
}