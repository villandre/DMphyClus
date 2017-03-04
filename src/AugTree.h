#include "TreeNode.h"

using namespace arma ;
using namespace Rcpp ;

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
};