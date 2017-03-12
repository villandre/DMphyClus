#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  
  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChild(TreeNode *) {assert(false) ;}
  std::vector<TreeNode *> GetChildren() {std::vector<TreeNode *> myVec; myVec.push_back(NULL) ; return myVec;} // An input node returns a null pointer when it is asked to provide the address of a child.
  void RemoveChildren() {}
  
  bool IsSolved() {return true ;}
  bool CanSolve() {return true ;}
  void SetSolved(bool status) {}
  void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, const uint &) {assert(false) ;}
  void ComputeSolution(solutionDictionaryType & dictionary, const mat &, const uint &, const uint &, const uint &) {assert(false) ;} //Solution is known, this should not get called.
  void InvalidateSolution() {assert(false) ;}
  vec GetSolution(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(std::floor(elementNum/numRateCats))->second.first ;};
  
  void SetInput(std::vector<uvec> * inputVec) { _inputVec = inputVec ;}
  std::vector<uvec> * GetInput() { return _inputVec ;}
  
  void InitMapAndIterVec(solutionDictionaryType &) ;
  void CopyIterVecAndExp() {assert(false) ; }
  void RestoreIterVecAndExp() {} // Solutions for input nodes are trivial and never change. It follows that a restore should not do anything.
  //std::vector<bool> UpdateDictionaryIter(solutionDictionaryType &, uint &) {assert(false) ; return std::vector<bool>(0) ;};
  mapIterator GetDictionaryIterator(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum) ;}
  S GetSfromVertex(const uint &, const uint &, const uint &) {assert(false) ;};
  
  fvec GetExponentIncrementVec(const uint & numRateCats) {return fvec(_dictionaryIterVec.size()*numRateCats, fill::zeros) ;}
  
  InputNode(uint & numLoci, uint & numRates)
  {
    _parent = NULL ;
    _dictionaryIterVec.reserve(numLoci*numRates) ;
    _updateFlag = false ;
  }
  
protected:
  std::vector<uvec> * _inputVec ;
};
