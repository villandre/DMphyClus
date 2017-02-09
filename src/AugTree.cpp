#include "AugTree.h"
#include<vector.h>
#include<memory.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(IntegerMatrix & edgeMatrix, List & alignmentBin, List & transProbMatList, NumericVector & limProbs) 
{
  _limProbs = as<vec>(limProbs) ;
  _transProbMatrixVec = as<std::vector<mat>>(transProbMatList) ;
  _alignmentBin = as<std::vector<longVec>>(alignmentBin) ;
  umat edgeMatrixCopy(as<umat>(edgeMatrix)) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
}

void AugTree::BuildTree(umat & edgeMatrix) 
{
  uint numTips = _alignmentBin[0].size();
  _tree.resize(edgeMatrix.n_rows + 1) ;
  for (umat::iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end(); iter = iter + 2) { // Iterator will move beyond the end of the matrix...
    if (*(iter+1) < numTips) 
    {
      InputNode myNode = InputNode() ;
      _tree.at(*iter) = &myNode; 
    } 
    else 
    {
      IntermediateNode myNode = IntermediateNode() ;
      _tree.at(*iter) = &myNode ; 
    }
    _tree.at(*iter)->SetId(*(iter + 1)) ;
  }
  // The vertices are all disjoint, the next loop defines their relationships
  for (umat::iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end(); iter = iter + 2) 
  { 
    _tree.at(*iter)->AddChild(_tree.at(*(iter+1))) ;
    _tree.at(*(iter+1))->SetParent(_tree.at(*iter)) ;
  }
}

SEXP AugTree::BuildEdgeMatrix() 
{
  umat edgeMatrix(_tree.size(), 2) ;
  std::transform(_tree.begin(), _tree.end(), edgeMatrix.begin(), [] (TreeNode * NodePointer) { return NodePointer->GetParent() ;}) ;
  return wrap(edgeMatrix) ;
}

long double AugTree::ComputeLogLik() 
{
  for (auto & i : _alignmentBin) 
  {
    InitializeTips(i) ;
    for (auto tup : zip(_withinTransProbMatrixVec, _betweenTransProbMatrix))
    SolveRoot() ;
  }
  
}

void AugTree::InitializePatterns(TreeNode * startNode) 
{
  for (auto & i : startNode->_children)
  {
    if (i->CanSolve())
    {
      i->SetPattern() ;
      i->ToggleSolved() ;
    }
    else
    {
     InitializePatterns(i) ;
    }
  }
  // The nodes had been marked as solved only so that we could get a pattern indicator for each. They are not really solved and so, must be marked as such.
  // Note that, as it should, ToggleSolved() does nothing for InputNodes.
  for (auto & i : _tree)
  {
    i->ToggleSolved() ;
  }
}

bool AugTree::SolveRoot() {
  TreeNode * currentPoint = _tree.at(_alignmentBin[0].size()) ;
  if (_dictionary.find(currentPoint->GetPattern()) == _dictionary.end()) {
    
  }
}