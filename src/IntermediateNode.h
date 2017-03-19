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
  void ComputeSolutions(solutionDictionaryType & solutionDictionary, const std::vector<mat> & transProbMats, const uint & transMatIndex, ThreadPool * myThreadpool) ;
  void ComputeSolution(mapIterator & solutionIter, solutionDictionaryType & solutionDictionary, const std::vector<mat> & transProbMatVec, const uint & locusNum, const uint & transMatrixIndex) ;
  void InvalidateSolution() ;
  
  void SetInput(std::vector<uvec> *) { assert(false) ;};
  std::vector<uvec> * GetInput() {assert(false) ; return NULL;} ;
  
  void InitMapAndIterVec(solutionDictionaryType &, const uint &) {assert(false) ;};
  
  void RestoreIterVecAndExp() 
  {
    if (_updateFlag) 
    {
      std::copy(_previousIterVec.begin(), _previousIterVec.end(), _dictionaryIterVec.begin()) ;
    }
    _updateFlag = false ;
  }
  mapIterator GetDictionaryIterator(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum) ;}
  S GetSfromVertex(const uint &, const uint &, const uint &) ;
  
  void PrepareSchedule(const solutionDictionaryType & solutionDictionary, iterVec & iteratorVec, const uint & locusNum, const uint & transMatrixIndex, const uint & numRates) ;
  
  IntermediateNode(uint & numLoci, uint & numRates, solutionDictionaryType & solutionDictionary): _isSolved(false) {
    _parent = NULL ;
    _dictionaryIterVec.resize(numLoci) ;
    _previousIterVec.resize(numLoci) ;
    
    for (uint i = 0 ; i < numLoci ; i++)
    {
      _dictionaryIterVec.at(i) = solutionDictionary->begin() ;
    }
    _updateFlag = false ;
  };
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   iterVec _previousIterVec ;
};

#endif /* INTERMEDIATENODE_H */
