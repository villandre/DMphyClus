#include "TreeNode.h"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

using namespace arma ;
using namespace Rcpp ; // Tandy Warnow

template <typename T>
void deallocate_container(T& c)
{
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i; 
}

class AugTree
{
protected:
  std::vector<TreeNode *> _tree ;
  mat _withinTransProbMatrix ;
  mat _betweenTransProbMatrix ;
  Col<double> _limProbs ;
  double _likelihood ;
  uint _numTips ;
  uint _rateCateg ;

  void BuildTree(umat &) ;
  void SolveOneLevel() ;
  void InitializeFromDictionary() ;
  void InitializeTips(const std::vector<uvec> &) ;
  void AssociateTransProbMatrices(const uvec &, const mat &, const mat &) ;
  void BindMatrixChildren(TreeNode *, const mat &, const bool) ;
  void PatternLookup(solutionDictionaryType &, TreeNode *) ;

public:
  AugTree(const umat &, const uvec &, const mat &, const mat &, const std::vector<uvec> &, const Col<double> &, const uint, const uint, solutionDictionaryType &) ;
  void TrySolve(TreeNode *)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &) ;
  SEXP BuildEdgeMatrix() ;
  void ComputeKeys(TreeNode *, solutionDictionaryType &) ;
  double GetLikelihood() const {return _likelihood ;} ;
  mat GetWithinTransProbMatrix() const {return _withinTransProbMatrix ;} ;
  mat GetBetweenTransProbMatrix() const {return _betweenTransProbMatrix ;} ;
  void SetWithinTransProbMatrix(mat withinTransProbs) {_withinTransProbMatrix = withinTransProbs;} ;
  void SetBetweenTransProbMatrix(mat betweenTransProbs) {_betweenTransProbMatrix = betweenTransProbs ;} ;
  ~AugTree() {deallocate_container(_tree) ;};
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
  Forest(const IntegerMatrix &, const NumericVector &, const List &, const List &, const List &, const NumericVector &, const uint, solutionDictionaryType &) ;
  void ComputeLoglik() ;
  double GetLoglik() {return _loglik ;} ;
  void NNmovePropagate() ;
};