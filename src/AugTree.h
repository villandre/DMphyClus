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
  mat _withinTransProbMatrix ;
  mat _betweenTransProbMatrix ;
  Col<long double> _limProbs ;
  double _likelihood ;
  uint _numTips ;
  
  void BuildTree(umat &) ;
  void SolveOneLevel() ;
  void SetPatterns() ;
  void InitializeFromDictionary() ;
  void InitializeTips(const std::vector<uvec> &) ;
  void AssociateTransProbMatrices(const NumericVector &, const mat &, const mat &, const uint) ;
  void BindMatrixChildren(TreeNode *, const mat &, const uint, const bool) ;
  void PatternLookup(solutionDictionaryType &, TreeNode *) ;
  
public:
  AugTree(const IntegerMatrix &, const NumericVector &, const mat &, const mat &, const std::vector<uvec> &, const NumericVector &, const uint, const uint) ;
  void TrySolve(TreeNode *)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &) ;
  SEXP BuildEdgeMatrix() ;
  void IdentifyPatterns(TreeNode *) ;
  double GetLikelihood() const {return _likelihood ;} ;
  mat GetWithinTransProbMatrix() const {return _withinTransProbMatrix ;} ;
  mat GetBetweenTransProbMatrix() const {return _betweenTransProbMatrix ;} ;
  void SetWithinTransProbMatrix(mat withinTransProbs) {_withinTransProbMatrix = withinTransProbs;} ;
  void SetBetweenTransProbMatrix(mat betweenTransProbs) {_betweenTransProbMatrix = betweenTransProbs ;} ;
};

class Forest 
{
protected:
  std::vector<AugTree> _forest ;
  double _loglik ;
  nodePatternDictionaryType _nodePatternDictionary ;
  solutionDictionaryType _solutionDictionary ; 
  uint _numLoci ;
  uint _numRateCats ; 

public:
  Forest(const IntegerMatrix &, const NumericVector &, const List &, const List &, const List &, const NumericVector &, const uint) ;
  void ComputeLoglik() ;
  double GetLoglik() {return _loglik ;} ;
  void NNmovePropagate() ;
};