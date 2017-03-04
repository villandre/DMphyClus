#ifndef FOREST_H
#define FOREST_H

#include "TreeNode.h"
#include "AugTree.h"

using namespace arma ;
using namespace Rcpp ; 

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

#endif /* FOREST_H */