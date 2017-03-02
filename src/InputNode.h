#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  bool IsSolved() {return true ;} ;
  bool CanSolve() {return true ;} ;
  bool CanFindKey() {return true ;};
  void AddChild(TreeNode * child) {assert(false) ;};
  void RemoveChildren() {};
  void RemoveChild(TreeNode *) {assert(false) ;} ;
  void SetSolution(vec & inputVec) {assert(false) ;};
  void ComputeSolution(solutionDictionaryType & dictionary, const mat &, double *) {assert(false) ;}; //Solution is known, this should not get called.
  void InvalidateSolution() {assert(false) ;};
  void SetSolved(bool status) {};
  void SetInput(uvec * inputVec) { _input = inputVec ;} ;
  std::vector<TreeNode *> GetChildren() {std::vector<TreeNode *> myVec; myVec.push_back(NULL) ; return myVec;}; // An input node returns a null pointer when it is asked to provide the address of a child.
  void DeriveKey(solutionDictionaryType &) ;
  vec GetSolution() {return conv_to<vec>::from(*_input) ;} ;
  uvec * GetInput() { return _input ;} ;
  void EnterSolution(TreeNode * originVertex) 
  {
    _input = originVertex->GetInput() ;
  };
  
  
  
  InputNode() {_parent = NULL ; _keyDefined = false ;};
  

protected:
  uvec * _input ;
};
