#include "AugTree.h"

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(IntegerMatrix & edgeMatrix, IntegerMatrix & alignmentBin, List & transProbMatList, NumericVector & limProbs) 
{
  _limProbs = as<vec>(limProbs) ;
  _transProbMatrixVec = as<std::vector<mat>>(transProbMatList) ;
  _alignmentBin = as<umat>(alignmentBin) ;
  umat edgeMatrixCopy(as<umat>(edgeMatrix)) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
}

void AugTree::BuildTree(umat & edgeMatrix) 
{
  std::vector<InputNode>() ;
  
}