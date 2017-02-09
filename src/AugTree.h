#include "TreeNode.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include <RcppArmadillo.h>

typedef std::vector<Col<long double>> longVec ;

using namespace arma ;
using namespace Rcpp ; // Tandy Warnow

class AugTree 
{
protected:
  std::vector<TreeNode *> _tree ;
  std::vector<mat> _withinTransProbMatrixVec ;
  std::vector<mat> _betweenTransProbMatrixVec ;
  vec _limProbs ;
  std::vector<longVec> _alignmentBin ;
  std::unordered_map<mapKey, Col<long double>, MyHash> _dictionary ;
  
  void BuildTree(umat &) ;
  
public:
  AugTree(IntegerMatrix &, List &, List &, NumericVector &) ;
  bool TrySolve(TreeNode *)  ;
  long double ComputeLogLik() ;
  void NearestNeighbourSwap() ;
  bool SolveRoot() ;
  SEXP BuildEdgeMatrix() ;
  void InitializePatterns(TreeNode *) ;
  void InitializeTips(longVec &) ;
};