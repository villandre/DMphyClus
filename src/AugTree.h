#include "TreeNode.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include <RcppArmadillo.h>

using namespace arma ;
using namespace Rcpp ; // Tandy Warnow

typedef std::unordered_map<std::size_t, Col<long double>, MyHash> dictionary ;

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
  void AssociateTransProbMatrices(const NumericVector &, const mat &, const mat &) ;
  void BindMatrixChildren(TreeNode *, const mat &) ;
  void PatternLookup(dictionary &, TreeNode *) ;
  
public:
  AugTree(const IntegerMatrix &, const NumericVector &, const mat &, const mat &, const std::vector<uvec> &, const NumericVector &, const uint) ;
  void TrySolve(TreeNode *)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(dictionary &) ;
  SEXP BuildEdgeMatrix() ;
  void IdentifyPatterns(TreeNode *) ;
  double GetLikelihood() const {return _likelihood ;} ;
  mat GetWithinTransProbMatrix() {return _withinTransProbMatrix ;} ;
  mat GetBetweenTransProbMatrix() {return _betweenTransProbMatrix ;} ;
  void SetWithinTransProbMatrix(mat withinTransProbs) {_withinTransProbMatrix = withinTransProbs;} ;
  void SetBetweenTransProbMatrix(mat betweenTransProbs) {_betweenTransProbMatrix = betweenTransProbs ;} ;
};

class Forest 
{
protected:
  std::vector<AugTree> _forest ;
  double _loglik ;
  dictionary _dictionary ;
  uint _numLoci ;
  uint _numRateCats ; 

public:
  Forest(const IntegerMatrix &, const NumericVector &, const List &, const List &, const List &, const NumericVector &, const uint) ;
  void ComputeLoglik() ;
  double GetLoglik() {return _loglik ;} ;
};