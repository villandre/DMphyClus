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
  void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, const uint &, const std::vector<bool> &) ;
  void ComputeSolution(solutionDictionaryType &, const mat &, const uint &, const uint &, const uint &) ;
  void InvalidateSolution() ;
  vec GetSolution(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum)->second ;} ;
  
  void SetInput(std::vector<uvec> *) { assert(false) ;};
  std::vector<uvec> * GetInput() {assert(false) ; return NULL;} ;
  
  void InitMapAndIterVec(solutionDictionaryType &) {assert(false) ;};
  void CopyIterVecAndExp() 
  { 
    std::copy(_dictionaryIterVec.begin(), _dictionaryIterVec.end(), _previousIterVec.begin()) ;
    _previousExponentIncrementVec = _exponentIncrementVec ;
    _updateFlag = true ;
  }
  void RestoreIterVecAndExp() 
  {
    if (_updateFlag) 
    {
      std::copy(_previousIterVec.begin(), _previousIterVec.end(), _dictionaryIterVec.begin()) ;
      _exponentIncrementVec = _previousExponentIncrementVec ;
    }
    _updateFlag = false ;
  }
  std::vector<bool> UpdateDictionaryIter(solutionDictionaryType &, uint &) ;
  mapIterator GetDictionaryIterator(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum) ;}
  S GetSfromVertex(const uint &, const uint &, const uint &) ;
  
  vec GetExponentIncrementVec(const uint & numRateCats) {return _exponentIncrementVec ;}
  
  IntermediateNode(uint & numLoci, uint & numRates): _isSolved(false) {
    _parent = NULL ;
    int numElements = numLoci*numRates ;
    _dictionaryIterVec.resize(numElements) ;
    _previousIterVec.resize(numElements) ;
    _updateFlag = false ;
    _exponentIncrementVec = vec(numElements, fill::zeros) ;
    _previousExponentIncrementVec = vec(numElements, fill::zeros) ;
  };
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   vec _exponentIncrementVec ;
   vec _previousExponentIncrementVec ;
};

#endif /* INTERMEDIATENODE_H */
