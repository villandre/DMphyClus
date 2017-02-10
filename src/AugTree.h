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
  std::vector<mat> _withinTransProbMatrixVec ;
  std::vector<mat> _betweenTransProbMatrixVec ;
  vec _limProbs ;
  std::vector<longVec> _alignmentBin ;
  std::unordered_map<mapKey, Col<long double>, MyHash> _dictionary ;
  double _logLik ;
  uint _numRateCats ;
  
  void BuildTree(umat &) ;
  void SolveOneLevel(uint locusNum, uint rateNum) ;
  void SetPatterns() ;
  void InitializeFromDictionary(uint, uint) ; 
  
public:
  AugTree(IntegerMatrix &, List &, List &, List &, NumericVector &) ;
  bool TrySolve(TreeNode *)  ;
  void NearestNeighbourSwap() ;
  bool SolveRoot() ;
  SEXP BuildEdgeMatrix() ;
  void InitializePatterns(TreeNode *) ;
  void InitializeTips(longVec &) ;
  double GetLogLik() {return _logLik ;} ;
  std::vector<mat> GetWithinTransProbMatrix() {return _withinTransProbMatrixVec ;} ;
  std::vector<mat> GetBetweenTransProbMatrix() {return _betweenTransProbMatrixVec ;} ;
  std::vector<longVec> GetLogSolutions() ; // R cannot understand long doubles. We'll have to work on the log scale.
};

class Forest 
{
private:
  std::vector<AugTree> _forest ;

public:
  Forest(std::vector<longVec> &, )
  
};