#include "TreeNode.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include <RcppArmadillo.h>

using namespace arma ;
using namespace Rcpp ;

class AugTree 
{
protected:
  std::vector<TreeNode *> _tree ;
  std::vector<mat> _transProbMatrixVec ;
  vec _limProbs ;
  std::unordered_map<int, Col<long double>> dictionary ;
  umat _alignmentBin ;
  
public:
  AugTree(IntegerMatrix &, IntegerMatrix &, List &, NumericVector &) ;
  void BuildTree(umat &) ;
  void ComputeLik() ;
};