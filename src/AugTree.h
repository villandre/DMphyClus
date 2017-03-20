#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include "TreeNode.h"
#include "ThreadPoolOther.h"

using namespace arma ;
using namespace Rcpp ;

class AugTree
{
protected:
  boost::shared_mutex _mutex ;
  unsigned int _numThreads ;
  ThreadPool * _threadpool ;
  
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
  
  // vec _likPropVec ; // This is scaled to avoid computational zeros.

  void BuildTree(const umat &) ;
  void SolveOneLevel() ;
  void InitializeFromDictionary() ;
  void InitializeVertices() ;
  void AssociateTransProbMatrices(const uvec &) ;
  void GetNNIverticesInternalWithin(TreeNode *, std::vector<uint> *) ;
  void GetNNIverticesInternalBetween(TreeNode *, std::vector<uint> *, uvec &) ;
  void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  
public:
  AugTree(const umat & edgeMatrix, const uvec & clusterMRCAs, std::vector<std::vector<uvec>> * alignmentBin, solutionDictionaryType & solutionDictionary, const uint & withinMatListIndex, const uint & betweenMatListIndex, const uint & numRateCats, gsl_rng * RNGpoint, unsigned int & numThreads) ;
  void BuildTreeNoAssign(const umat &) ;
  void TrySolve(TreeNode *, const std::vector<mat> &, const std::vector<mat> &, const uint &)  ;
  void NearestNeighbourSwap() ;
  void SolveRoot(solutionDictionaryType &, const mat &, const mat &, const vec &, const uint &) ;
  
  void BindMatrixBetween(TreeNode *, const mat &) ;
  void InvalidateAll() ;
  void BindMatrix(TreeNode *, const bool) ;
  umat BuildEdgeMatrix() ;
  std::vector<uint> GetTwoVerticesForNNI(TreeNode *, uvec &) ;
  void CopyAugTreeNonPointer(AugTree *) ;
  void CheckAndInvalidateBetweenRecursive(TreeNode *) ; 
  void RestorePreviousConfig(const IntegerMatrix &, const bool, const int &, const int &, const IntegerVector &, bool) ;
  void NegateAllUpdateFlags() ;
  
  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;
  std::vector<uint> GetNNIverticesWithin(TreeNode *) ;
  std::vector<uint> GetNNIverticesBetween(TreeNode *, uvec &) ;

  void ComputeLoglik(const std::vector<mat> &, const std::vector<mat> &, const vec &) ;
  double GetLoglik() {return _logLik ;}
  gsl_rng * GetRandomNumGenerator() {return _randomNumGenerator ;}
  
  uint GetNumRateCats() {return _numRateCats ;}
  uint GetNumLoci() {return _numLoci ;}
  std::vector<std::vector<uvec>> * GetAlignmentBinReference() {return _alignmentBinReference ;}
  uint GetNumTips() {return _numTips ;}
  solutionDictionaryType GetSolutionDictionary() { return _solutionDictionary ;}
  uint GetWithinMatListIndex() {return _withinMatListIndex ;}
  uint GetBetweenMatListIndex() {return _betweenMatListIndex ;}
  unsigned int GetNumThreads() {return _numThreads ;}
  
  void InvalidateBetweenSolutions() ;
  void InvalidateAllSolutions() ;
 
  void SetBetweenMatListIndex(const uint & index) {_betweenMatListIndex = index ;}
  void SetWithinMatListIndex(const uint & index) {_withinMatListIndex = index ;}
  void HandleSplit(uint) ;
  void HandleMerge(uvec &) ;
  void SetLogLik(double logLik) {_logLik = logLik ;}
  void RebuildTrees(const umat &) ;
  void SetRNG(gsl_rng * myRNG) { _randomNumGenerator = myRNG ;}
  
  void RearrangeTreeNNI(uint, uint) ;
  
  //void ComputeLoglik(List &, List &, NumericVector &) ;
  void PrintSolutions(const uint &) ;
  
  void DestroyThreads()
  {
    _threadpool->JoinAll() ;
    delete _threadpool ;
  }
  
  ~AugTree() {deallocate_container(_vertexVector) ;};
};