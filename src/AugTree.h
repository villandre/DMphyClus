#include "TreeNode.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include <RcppArmadillo.h>

using namespace arma ;
using namespace Rcpp ; // Tandy Warnow

class AugTree 
{
protected:
  std::vector<TreeNode *> _tree ;
  std::vector<mat> _transProbMatrixVec ;
  vec _limProbs ;
  std::unordered_map<int, Col<long double>> dictionary ;
  std::vector<umat> _alignmentBin ;
  
  void BuildTree(umat &) ;
  
public:
  AugTree(IntegerMatrix &, List &, List &, NumericVector &) ;
  bool TrySolve(TreeNode *)  ;
  long double ComputeLogLik() ;
  void NearestNeighbourSwap() ;
  bool SolveRoot() ;
  SEXP BuildEdgeMatrix() ;
};