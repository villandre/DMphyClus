#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  bool IsSolved() {return true ;} ;
  bool CanSolve() {return true ;} ;
  void AddChild(TreeNode* child) {assert(false) ;};
  void RemoveChild(TreeNode* child) {assert(false) ;};
  void SetSolution(Col<long double> & inputVec) {assert(false) ;};
  void GetSolution() const;
  void ComputeSolution() {assert(false) ;}; //Solution is known, this should not get called.
  void ToggleSolved() {};
  void SetPattern() ; 
  void SetInput(const Col<bool> & inputVec) { _input = inputVec ;} ;
  void DeriveKey(solutionDictionaryType &) ;
  Col<long double> GetSolution() {return conv_to<Col<long double>>::from(_input) ;} ;
  std::size_t GetPattern() {assert(false); return 0 ;}; // Should not get called, since tips are by default resolved at the start.
  std::vector<TreeNode *> GetChildren() {return std::vector<TreeNode *>{NULL} ;}; // An input node returns a null pointer when it is asked to provide the address of a child.
  
  InputNode() {}; 
  
protected:
  void InvalidateSolution() {assert(false) ;};
  Col<bool> _input ;
};