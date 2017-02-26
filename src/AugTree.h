#include "TreeNode.h"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#include <gsl/gsl_rng.h>

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
  std::vector<TreeNode *> _vertexVector ;
  mat _withinTransProbMatrix ;
  mat _betweenTransProbMatrix ;
  Col<double> _limProbs ;
  double _likelihood ;
  uint _numTips ;
  uint _rateCateg ;

  void BuildTree(umat &) ;
  void SolveOneLevel() ;
  void InitializeFromDictionary() ;
  void InitializeVertices(const std::vector<uvec> &) ;
  void AssociateTransProbMatrices(const uvec &, const mat &, const mat &) ;
  void PatternLookup(solutionDictionaryType &, TreeNode *) ;
  void GetNNIverticesInternalWithin(TreeNode *, std::vector<uint> *) ;
  void GetNNIverticesInternalBetween(TreeNode *, std::vector<uint> *, uvec &) ;
  void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  
public:
  AugTree(const umat &, const uvec &, const mat &, const mat &, const std::vector<uvec> &, const Col<double> &, const uint, const uint, solutionDictionaryType &) ;
  AugTree(const umat &, const vec &, const uint, const uint) ;
  
  void TrySolve(TreeNode *, solutionDictionaryType &)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &) ;
  void ComputeKeys(TreeNode *, solutionDictionaryType &) ;
  void BindMatrixBetween(TreeNode *, const mat &) ;
  void InvalidateAll() ;
  void BindMatrix(TreeNode *, const mat &, const bool) ;
  umat BuildEdgeMatrix() ;
  std::vector<uint> GetTwoVerticesForNNI(gsl_rng *, TreeNode *, uvec &) ;
  void CopyAugTreeNonPointer(AugTree *) ;
  
  void SetWithinTransProbMatrix(mat withinTransProbs) {_withinTransProbMatrix = withinTransProbs ;} ;
  void SetBetweenTransProbMatrix(mat betweenTransProbs) {_betweenTransProbMatrix = betweenTransProbs ;} ;
  
  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;
  mat GetWithinTransProbMatrix() const {return _withinTransProbMatrix ;} ;
  mat GetBetweenTransProbMatrix() const {return _betweenTransProbMatrix ;} ;
  double GetLikelihood() const {return _likelihood ;} ;
  uint GetNumTips() {return _numTips ;} ;
  std::vector<uint> GetNNIverticesWithin(TreeNode *) ;
  std::vector<uint> GetNNIverticesBetween(TreeNode *, uvec &) ;
  vec GetLimProbs() { return _limProbs ;} ;
  
  void RearrangeTreeNNI(uint, uint) ;
  
  ~AugTree() {deallocate_container(_vertexVector) ;};
  // AugTree( const AugTree& other ):_limProbs(other._limProbs), _likelihood(other._likelihood), _numTips(other._numTips), _rateCateg(other._rateCateg)
  // {
  //   _withinTransProbMatrix = other._withinTransProbMatrix ;
  //   _betweenTransProbMatrix = other._betweenTransProbMatrix ;
  //   std::transform(other._vertexVector.begin(), other._vertexVector.end(), std::back_inserter(_vertexVector), std::mem_fn(&TreeNode::clone)) ;
  // }; // A copy constructor... The default copy constructor will copy the addresses in _vertexVector, rather than allocating a new vector and copying the elements pointed to by _vertexVector...
};

class Forest
{
protected:
  std::vector<AugTree *> _forest ;
  double _loglik ;
  //nodePatternDictionaryType _nodePatternDictionary ;
  solutionDictionaryType _solutionDictionary ;
  uint _numLoci ;
  uint _numRateCats ;
  gsl_rng * _randomNumGenerator ;

public:
  Forest(const IntegerMatrix &, const NumericVector &, const List &, const List &, const List &, const NumericVector &, const uint, const uint, solutionDictionaryType &) ;
  Forest() ;
  Forest(const IntegerMatrix &, const vec &, uint, uint, uint, gsl_rng *, solutionDictionaryType) ;
  
  ~Forest() {deallocate_container(_forest) ;}
  // Forest( const Forest& other ):_loglik(other._loglik), _numLoci(other._numLoci), _numRateCats(other._numRateCats)
  // {
  //   std::transform(other._forest.begin(), other._forest.end(), std::back_inserter(_forest), std::mem_fn(&AugTree::clone)) ;
  // }; // A copy constructor... The default copy constructor will copy the addresses in _vertexVector, rather than allocating a new vector and copying the elements pointed to by _vertexVector...
  // 
  // I did not explicitly specify a copy constructor, because the default one will call the copy constructor on all members, including the one I defined for AugTree, which produces deep copies of _vertexVector.
  //void ForestDeepCopy(Forest *) ;
  
  void ComputeLoglik() ;
  double GetLoglik() {return _loglik ;} ;
  gsl_rng * GetRandomNumGenerator() {return _randomNumGenerator ;} ;
  std::vector<AugTree *> GetForest() {return _forest ;} ;
  uint GetNumRateCats() {return _numRateCats ;} ;
  uint GetNumLoci() {return _numLoci ;}
  solutionDictionaryType GetSolutionDictionary() { return _solutionDictionary ;} ;
  
  void AmendBetweenTransProbs(std::vector<mat> &) ;
  void AmendWithinTransProbs(std::vector<mat> &, uvec &) ;
  void HandleSplit(uint, std::vector<mat> &) ;
  void HandleMerge(uvec &, std::vector<mat> &) ;
  void SetLogLik(double logLik) {_loglik = logLik ;} ;
  
  void InputForestElements(Forest * originForest) ;
};