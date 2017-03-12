#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H

class IntermediateNode : public TreeNode
{
public:
  void AddChild(TreeNode * child) {_children.push_back(child) ;};
  void RemoveChild(TreeNode *) ;
  std::vector<TreeNode *> GetChildren() {return _children;};
  void RemoveChildren() {_children.clear() ;} ;
  
  bool IsSolved() {return _isSolved ;};
  bool CanSolve() ;
  void SetSolved(bool status) {_isSolved = status ;};
  void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, const uint &) ;
  void ComputeSolution(solutionDictionaryType &, const mat &, const uint &, const uint &, const uint &) ;
  void InvalidateSolution() ;
  vec GetSolution(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum)->second.first ;} ;
  
  void SetInput(std::vector<uvec> *) { assert(false) ;};
  std::vector<uvec> * GetInput() {assert(false) ; return NULL;} ;
  
  void InitMapAndIterVec(solutionDictionaryType &) {assert(false) ;};
  void CopyIterVecAndExp() 
  { 
    //std::copy(_dictionaryIterVec.begin(), _dictionaryIterVec.end(), _previousIterVec.begin()) ;
    _previousExponentIncrementVec = _exponentIncrementVec ;
  }
  void RestoreIterVecAndExp() 
  {
    if (_updateFlag) 
    {
      for (uint i = 0 ; 0 < _iterMove.size(); i++)
      {
        std::advance(_dictionaryIterVec.at(i), _iterMove.at(i)) ;
      }
      //std::fill(_iterMove.begin(), _iterMove.end(), 0) ;
      //std::copy(_previousIterVec.begin(), _previousIterVec.end(), _dictionaryIterVec.begin()) ;
      _exponentIncrementVec = _previousExponentIncrementVec ;
    }
    _updateFlag = false ;
  }
  //std::vector<bool> UpdateDictionaryIter(solutionDictionaryType &, uint &) ;
  mapIterator GetDictionaryIterator(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum) ;}
  S GetSfromVertex(const uint &, const uint &, const uint &) ;
  
  fvec GetExponentIncrementVec(const uint & numRateCats) {return _exponentIncrementVec ;}
  
  IntermediateNode(uint & numLoci, uint & numRates, solutionDictionaryType & solutionDictionary): _isSolved(false) {
    _parent = NULL ;
    int numElements = numLoci*numRates ;
    _dictionaryIterVec.resize(numElements) ;
    uint rateCateg = 0 ;
    for (uint i = 0 ; i < numElements ; i++)
    {
      _dictionaryIterVec.at(i) = solutionDictionary->at(rateCateg).begin() ;
      rateCateg = littleCycle(rateCateg+1, solutionDictionary->size()) ;
    }
    _iterMove.resize(numElements) ;
    std::fill(_iterMove.begin(), _iterMove.end(), 0) ;
    _updateFlag = false ;
    _exponentIncrementVec = fvec(numElements, fill::zeros) ;
    _previousExponentIncrementVec = fvec(numElements, fill::zeros) ;
  };
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   fvec _exponentIncrementVec ;
   fvec _previousExponentIncrementVec ;
};

#endif /* INTERMEDIATENODE_H */
