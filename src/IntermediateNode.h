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
  void ComputeSolution(solutionDictionaryType &, const std::vector<mat> &, const uint &, const uint &) ;
  void InvalidateSolution() ;
  
  void SetInput(std::vector<uvec> *) { assert(false) ;};
  std::vector<uvec> * GetInput() {assert(false) ; return NULL;} ;
  
  void InitMapAndIterVec(solutionDictionaryType &, const uint &) {assert(false) ;};
  
  void RestoreIterVecAndExp(solutionDictionaryType & solutionDictionary) 
  {
    if (_updateFlag) 
    {
      _dictionaryIterVec = _previousIterVec ;
    }
    _updateFlag = false ;
  }
  mapIterator GetDictionaryIterator(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum) ;}
  S GetSfromVertex(const uint &, const uint &, const uint &) ;
  
  IntermediateNode(uint & numLoci, uint & numRates, solutionDictionaryType & solutionDictionary): _isSolved(false) {
    _parent = NULL ;
    _dictionaryIterVec.resize(numLoci) ;
    for (uint i = 0 ; i < numLoci ; i++)
    {
      _dictionaryIterVec.at(i) = solutionDictionary->begin() ;
    }
    _updateFlag = false ;
  };
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
};

#endif /* INTERMEDIATENODE_H */
