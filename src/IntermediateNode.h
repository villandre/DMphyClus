#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H

class IntermediateNode : public TreeNode
{
public:

  bool IsSolved() {return _isSolved ;};
  bool CanSolve() ;
  bool CanFindKey() ;
  void AddChild(TreeNode * child) {_children.push_back(child) ;};
  void RemoveChildren() {_children.clear() ;} ;
  void RemoveChild(TreeNode *) ;
  void SetSolution(vec & solution) { _solution = solution ;};
  void ComputeSolution(solutionDictionaryType &, const mat &, double *) ;
  void InvalidateSolution() ;
  void SetSolved(bool status) {_isSolved = status ;};
  void SetInput(uvec *) { assert(false) ;};
  std::vector<TreeNode *> GetChildren() {return _children;};
  void DeriveKey(solutionDictionaryType &) ;
  vec GetSolution() {return _solution ;} ;
  uvec * GetInput() {assert(false) ; return NULL;} ;
  void EnterSolution(TreeNode * originVertex)
  {
    _solution = originVertex->GetSolution() ;
    _isSolved = originVertex->IsSolved() ;
  };

  IntermediateNode(): _isSolved(false) {_parent = NULL ;  _keyDefined = false ;};
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   Col<double> _solution ;
};

#endif /* INTERMEDIATENODE_H */
