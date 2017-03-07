#include "TreeNode.h"

using namespace arma ;
using namespace Rcpp ;

class AugTree
{
protected:
  double _logLik ;
  std::vector<TreeNode *> _vertexVector ;
  solutionDictionaryType _solutionDictionary ;
  uint _numLoci ;
  uint _numTips ;
  uint _numRateCats ;
  uint _withinMatListIndex ;
  uint _betweenMatListIndex ;
  gsl_rng * _randomNumGenerator ;
  
  std::vector<std::vector<uvec>> * _alignmentBinReference ;
  
  vec _likPropVec ; // This is scaled to avoid computational zeros.
  
  vec _exponentVec ;

  void BuildTree(const umat &) ;
  void SolveOneLevel() ;
  void InitializeFromDictionary() ;
  void InitializeVertices() ;
  void AssociateTransProbMatrices(const uvec &) ;
  void GetNNIverticesInternalWithin(TreeNode *, std::vector<uint> *) ;
  void GetNNIverticesInternalBetween(TreeNode *, std::vector<uint> *, uvec &) ;
  void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  
public:
  AugTree(const umat &, const uvec &, std::vector<std::vector<uvec>> *, solutionDictionaryType &, const uint &, const uint &, gsl_rng *) ;
  //AugTree(const umat &, const uint &, const uint &) ;
  
  void BuildTreeNoAssign(const umat &) ;
  void TrySolve(TreeNode *, solutionDictionaryType &, const mat &, const mat &)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &, const mat &, const mat &, const vec &, const uint &) ;
  void ComputeKeys(TreeNode *, const uint &, const uint &) ;
  void PatternLookup(solutionDictionaryType &, TreeNode *) ;
  void BindMatrixBetween(TreeNode *, const mat &) ;
  void InvalidateAll() ;
  void BindMatrix(TreeNode *, const bool) ;
  umat BuildEdgeMatrix(const uint &) ;
  std::vector<uint> GetTwoVerticesForNNI(gsl_rng *, TreeNode *, uvec &) ;
  void CopyAugTreeNonPointer(AugTree *) ;
  void CheckAndInvalidateBetweenRecursive(TreeNode *) ; 
  
  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;
  vec GetLikPropVec() const {return _likPropVec ;} ;
  std::vector<uint> GetNNIverticesWithin(TreeNode *) ;
  std::vector<uint> GetNNIverticesBetween(TreeNode *, uvec &) ;
  vec GetExponentVec() { return _exponentVec ;} ;
  void ComputeLoglik(List &, List &, const vec &) ;
  double GetLoglik() {return _logLik ;}
  gsl_rng * GetRandomNumGenerator() {return _randomNumGenerator ;}
  
  uint GetNumRateCats() {return _numRateCats ;}
  uint GetNumLoci() {return _numLoci ;}
  std::vector<std::vector<uvec>> * GetAlignmentBinReference() {return _alignmentBinReference ;}
  uint GetNumTips() {return _numTips ;}
  solutionDictionaryType GetSolutionDictionary() { return _solutionDictionary ;}
  uint GetWithinMatListIndex() {return _withinMatListIndex ;}
  uint GetBetweenMatListIndex() {return _betweenMatListIndex ;}
  void InvalidateBetweenSolutions() ;
  void InvalidateAllSolutions() ;
 
  void SetBetweenMatListIndex(const uint & index) {_betweenMatListIndex = index ;}
  void SetWithinMatListIndex(const uint & index) {_withinMatListIndex = index ;}
  void HandleSplit(uint) ;
  void HandleMerge(uvec &) ;
  void SetLogLik(double logLik) {_logLik = logLik ;}
  void RearrangeNNI(const uint, const uint) ;
  void RebuildTrees(const umat &) ;
  void SetRNG(gsl_rng * myRNG) { _randomNumGenerator = myRNG ;}
  
  void RearrangeTreeNNI(uint, uint, solutionDictionaryType) ;
  
  void ComputeLoglik(List &, List &, NumericVector &) ;
  
  ~AugTree() {deallocate_container(_vertexVector) ;};
};

class AugTreeWithUndoMove;