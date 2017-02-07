#include "AugTree.h"

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
  std::vector<InputNode>() ;
  
}