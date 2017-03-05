#include "TreeNode.h"

using namespace arma ;
using namespace Rcpp ;

class AugTree
{
protected:
  std::vector<TreeNode *> _vertexVector ;
  
  double _likelihoodProp ; // This is scaled to avoid computational zeros.
  uint _rateCateg ;
  double _exponentContainer ;

  void BuildTree(const umat &, const uint &) ;
  void SolveOneLevel() ;
  void InitializeFromDictionary() ;
  void InitializeVertices(std::vector<uvec> *, solutionDictionaryType &) ;
  void AssociateTransProbMatrices(const uvec &, const uint &) ;
  void GetNNIverticesInternalWithin(TreeNode *, std::vector<uint> *) ;
  void GetNNIverticesInternalBetween(TreeNode *, std::vector<uint> *, uvec &) ;
  void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  
public:
  AugTree(const umat &, const uvec &, std::vector<uvec> *, const uint, solutionDictionaryType &, const uint &) ;
  AugTree(const umat &, const uint &, const uint &) ;
  
  void BuildTreeNoAssign(const umat &) ;
  void TrySolve(TreeNode *, solutionDictionaryType &, const mat &, const mat &)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &, const mat &, const mat &, const vec &, const uint &) ;
  void ComputeKeys(TreeNode *, solutionDictionaryType &, const uint, const uint) ;
  void PatternLookup(solutionDictionaryType &, TreeNode *) ;
  void BindMatrixBetween(TreeNode *, const mat &) ;
  void InvalidateAll() ;
  void BindMatrix(TreeNode *, const bool) ;
  umat BuildEdgeMatrix(const uint &) ;
  std::vector<uint> GetTwoVerticesForNNI(gsl_rng *, TreeNode *, uvec &) ;
  void CopyAugTreeNonPointer(AugTree *) ;
  void CheckAndInvalidateBetweenRecursive(TreeNode *) ; 
  
  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;
  double GetLikelihood() const {return _likelihoodProp ;} ;
  std::vector<uint> GetNNIverticesWithin(TreeNode *) ;
  std::vector<uint> GetNNIverticesBetween(TreeNode *, uvec &) ;
  double GetExponentContainer() { return _exponentContainer ;} ;
  uint GetRateCateg() {return _rateCateg ;} ;
  
  void RearrangeTreeNNI(uint, uint, solutionDictionaryType) ;
  
  ~AugTree() {deallocate_container(_vertexVector) ;};
};