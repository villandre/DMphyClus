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
  double _likelihoodProp ; // This is scaled to avoid computational zeros.
  uint _numTips ;
  uint _rateCateg ;
  double _exponentContainer ;

  void BuildTree(umat &) ;
  void SolveOneLevel() ;
  void InitializeFromDictionary() ;
  void InitializeVertices(std::vector<uvec> *, solutionDictionaryType &) ;
  void AssociateTransProbMatrices(const uvec &) ;
  void GetNNIverticesInternalWithin(TreeNode *, std::vector<uint> *) ;
  void GetNNIverticesInternalBetween(TreeNode *, std::vector<uint> *, uvec &) ;
  void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  
public:
  AugTree(const umat &, const uvec &, std::vector<uvec> *, const Col<double> &, const uint, const uint, solutionDictionaryType &) ;
  AugTree(const umat &, const vec &, const uint, const uint) ;
  
  void TrySolve(TreeNode *, solutionDictionaryType &, const mat &, const mat &)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &, const mat &, const mat &) ;
  void ComputeKeys(TreeNode *, solutionDictionaryType &, const uint, const uint) ;
  void PatternLookup(solutionDictionaryType &, TreeNode *) ;
  void BindMatrixBetween(TreeNode *, const mat &) ;
  void InvalidateAll() ;
  void BindMatrix(TreeNode *, const bool) ;
  umat BuildEdgeMatrix() ;
  std::vector<uint> GetTwoVerticesForNNI(gsl_rng *, TreeNode *, uvec &) ;
  void CopyAugTreeNonPointer(AugTree *) ;
  void CheckAndInvalidateBetweenRecursive(TreeNode *) ; 
  
  void SetWithinTransProbMatrix(mat withinTransProbs) {_withinTransProbMatrix = withinTransProbs ;} ;
  void SetBetweenTransProbMatrix(mat betweenTransProbs) {_betweenTransProbMatrix = betweenTransProbs ;} ;
  
  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;
  mat GetWithinTransProbMatrix() const {return _withinTransProbMatrix ;} ;
  mat GetBetweenTransProbMatrix() const {return _betweenTransProbMatrix ;} ;
  double GetLikelihood() const {return _likelihoodProp ;} ;
  uint GetNumTips() {return _numTips ;} ;
  std::vector<uint> GetNNIverticesWithin(TreeNode *) ;
  std::vector<uint> GetNNIverticesBetween(TreeNode *, uvec &) ;
  vec GetLimProbs() { return _limProbs ;} ;
  double GetExponentContainer() { return _exponentContainer ;} ;
  uint GetRateCateg() {return _rateCateg ;} ;
  
  void RearrangeTreeNNI(uint, uint, solutionDictionaryType) ;
  
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
  uint _withinMatListIndex ;
  uint _betweenMatListIndex ;
  gsl_rng * _randomNumGenerator ;
  std::vector<mat> _withinTransProbMatVec ;
  std::vector<mat> _betweenTransProbMatVec ;
  std::vector<std::vector<uvec>> * _alignmentBinReference ;
  
public:
  Forest(const IntegerMatrix &, const NumericVector &, std::vector<std::vector<uvec>> *, const List &, const List &, const NumericVector &, const uint, const uint, solutionDictionaryType &, const uint, const uint) ;
  Forest() ;
  Forest(const IntegerMatrix &, const vec &, uint, uint, uint, gsl_rng *, solutionDictionaryType, std::vector<std::vector<uvec>> *, const uint, const uint) ;
  
  ~Forest() {deallocate_container(_forest) ;}
  
  void ComputeLoglik() ;
  double GetLoglik() {return _loglik ;} ;
  gsl_rng * GetRandomNumGenerator() {return _randomNumGenerator ;} ;
  std::vector<AugTree *> GetForest() {return _forest ;} ;
  uint GetNumRateCats() {return _numRateCats ;} ;
  uint GetNumLoci() {return _numLoci ;}
  std::vector<std::vector<uvec>> * GetAlignmentBinReference() {return _alignmentBinReference ;} ;
  
  solutionDictionaryType GetSolutionDictionary() { return _solutionDictionary ;} ;
  std::vector<mat> GetWithinTransProbMatVec() { return _withinTransProbMatVec ;} ;
  std::vector<mat> GetBetweenTransProbMatVec() { return _betweenTransProbMatVec ;} ;
  uint GetWithinMatListIndex() {return _withinMatListIndex ;} ;
  uint GetBetweenMatListIndex() {return _betweenMatListIndex ;} ;
  
  void InvalidateBetweenSolutions() ;
  void InvalidateAllSolutions() ;
  
  void SetBetweenTransProbs(const std::vector<mat> newBetweenTransProbsVec) {_betweenTransProbMatVec = newBetweenTransProbsVec ;} ;
  void SetWithinTransProbs(const std::vector<mat> newWithinTransProbsVec) {_withinTransProbMatVec = newWithinTransProbsVec ;} ;
  void HandleSplit(uint) ;
  void HandleMerge(uvec &) ;
  void SetLogLik(double logLik) {_loglik = logLik ;} ;
  void RearrangeNNI(const uint, const uint) ;
  
  void InputForestElements(XPtr<Forest> originForest) ;
};