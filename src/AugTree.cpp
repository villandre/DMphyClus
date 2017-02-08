#include "AugTree.h"
#include<vector.h>
#include<memory.h>

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(IntegerMatrix & edgeMatrix, List & alignmentBin, List & transProbMatList, NumericVector & limProbs) 
{
  _limProbs = as<vec>(limProbs) ;
  _transProbMatrixVec = as<std::vector<mat>>(transProbMatList) ;
  _alignmentBin = as<std::vector<umat>>(alignmentBin) ;
  umat edgeMatrixCopy(as<umat>(edgeMatrix)) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
}

void AugTree::BuildTree(umat & edgeMatrix) 
{
  uint numTips = _alignmentBin[0].n_cols ;
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
  for (auto &i : _alignmentBin)
  {
    cout << "hey" ;
  }
  return 0 ;
}