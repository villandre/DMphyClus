#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  bool IsSolved() {return true ;}
  bool CanSolve() {return true ;}
  bool CanFindKey() {return true ;}
  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChildren() {}
  void RemoveChild(TreeNode *) {assert(false) ;}
  void ComputeSolution(solutionDictionaryType & dictionary, const mat &, double &, const uint &, const uint &) {assert(false) ;} //Solution is known, this should not get called.
  void InvalidateSolution() {assert(false) ;}
  void SetSolved(bool status) {}
  void SetInput(std::vector<uvec> * inputVec) { _input = inputVec ;}
  std::vector<TreeNode *> GetChildren() {std::vector<TreeNode *> myVec; myVec.push_back(NULL) ; return myVec;} // An input node returns a null pointer when it is asked to provide the address of a child.
  void DeriveKey(solutionDictionaryType &, const uint &, const uint &, const uint &) ;
  vec GetSolution(solutionDictionaryType & dictionary, const uint & rateCateg, const std::size_t & dictionaryKey) {return conv_to<vec>::from(*_input) ;}
  std::vector<uvec> * GetInput() { return _input ;}
 // void EnterInput(TreeNode * originVertex) { _input = originVertex->GetInput() ;}
  void MarkKeyUndefined() {}
  
  InputNode(uint & numLoci)
  {
    _parent = NULL ; 
    _keyDefined = false ;
    _dictionaryKeyVec.reserve(numLoci) ;
  }
  
protected:
  std::vector<uvec> * _input ;
};
