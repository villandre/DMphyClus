#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  bool IsSolved() {return true ;}
  bool CanSolve() {return true ;}
  
  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChildren() {}
  void RemoveChild(TreeNode *) {assert(false) ;}
  void ComputeSolution(solutionDictionaryType & dictionary, const mat &, double &, const uint &, const uint &, const uint &) {assert(false) ;} //Solution is known, this should not get called.
  void InvalidateSolution() {assert(false) ;}
  void SetSolved(bool status) {}
  void SetInput(std::vector<uvec> * inputVec) { _inputVec = inputVec ;}
  std::vector<TreeNode *> GetChildren() {std::vector<TreeNode *> myVec; myVec.push_back(NULL) ; return myVec;} // An input node returns a null pointer when it is asked to provide the address of a child.
  
  vec GetSolution(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(std::floor(elementNum/3))->second ;};
  std::vector<uvec> * GetInput() { return _inputVec ;}
  void MarkKeyUndefined() {}
  void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, vec &, const uint &) {assert(false) ;}
  void InitMapAndIterVec(solutionDictionaryType &) ;
  void CopyIterVec() {assert(false) ; }
  void RestoreIterVec() {} // Solutions for input nodes are trivial and never change. It follows that a restore should not do anything.
  
  InputNode(uint & numLoci)
  {
    _parent = NULL ;
    _updateFlag = false ;
  }
  
protected:
  std::vector<uvec> * _inputVec ;
};
