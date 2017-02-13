#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H

class IntermediateNode : public TreeNode
{
public:

  bool IsSolved() {return _isSolved ;};
  bool CanSolve() ;
  void AddChild(TreeNode * child) {_children.push_back(child) ;};
  void RemoveChild(TreeNode *) ;
  void SetSolution(Col<double> & solution) { _solution = solution ;};
  void ComputeSolution() ;
  void InvalidateSolution() ;
  void ToggleSolved() {_isSolved = !_isSolved ;};
  void SetInput(const uvec &) { assert(false) ;};
  std::vector<TreeNode *> GetChildren() {return _children;};
  void DeriveKey(solutionDictionaryType &) ;
  Col<double> GetSolution() {return _solution ;} ;

  IntermediateNode(): _isSolved(false)  {_parent = NULL ;};

protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   Col<double> _solution ;
};

#endif /* INTERMEDIATENODE_H */