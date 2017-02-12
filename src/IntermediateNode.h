#include <assert.h>
#include "TreeNode.h"

class IntermediateNode : public TreeNode
{
public:
  
  bool IsSolved() {return _isSolved ;};
  bool CanSolve() ;
  std::size_t GetPattern() {return _pattern ;};
  void AddChild(TreeNode * child) {_children.push_back(child) ;};
  void RemoveChild(TreeNode*) ;
  void SetSolution(Col<long double> &) ;
  void ComputeSolution() ;
  void InvalidateSolution() ;
  //void SetPattern(std::unordered_map<mapKey, Col<long double>, MyHash> &) ;
  void ToggleSolved() {_isSolved = !_isSolved ;};
  void SetInput(const std::vector<uint> &) { assert(false) ;};
  void DeriveKey(solutionDictionaryType &) ;
  Col<long double> GetSolution() {return _solution ;} ;
  IntermediateNode(): _isSolved(false)  {_parent = NULL ;};
  std::vector<TreeNode *> GetChildren() {return _children;};
  
protected:
  
  bool _isSolved ;
  std::size_t _pattern ;
  std::vector<TreeNode *> _children ;
  Col<long double> _solution ;
};