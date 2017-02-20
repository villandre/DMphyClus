#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  bool IsSolved() {return true ;} ;
  bool CanSolve() {return true ;} ;
  bool CanFindKey() {return true ;};
  void AddChild(TreeNode * child) {assert(false) ;};
  void RemoveChildren() {};
  void SetSolution(vec & inputVec) {assert(false) ;};
  void ComputeSolution(solutionDictionaryType & dictionary) {assert(false) ;}; //Solution is known, this should not get called.
  void ToggleSolved() {};
  void SetInput(const uvec & inputVec) { _input = inputVec ;} ;
  std::vector<TreeNode *> GetChildren() {return std::vector<TreeNode *>{NULL} ;}; // An input node returns a null pointer when it is asked to provide the address of a child.
  void DeriveKey(solutionDictionaryType &) ;
  vec GetSolution() {return conv_to<vec>::from(_input) ;} ;

  InputNode() {_parent = NULL ; _keyDefined = false ;};

protected:
  void InvalidateSolution() {assert(false) ;};
  uvec _input ;
};
