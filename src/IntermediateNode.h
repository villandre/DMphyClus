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
  //void SetSolution(vec * solution) { _solution = solution ;};
  void ComputeSolution(solutionDictionaryType &, const mat &, double *, const uint &) ;
  void InvalidateSolution() ;
  void SetSolved(bool status) {_isSolved = status ;};
  void SetInput(uvec *) { assert(false) ;};
  std::vector<TreeNode *> GetChildren() {return _children;};
  void DeriveKey(solutionDictionaryType &, const uint &, const uint &) ;
  vec GetSolution(solutionDictionaryType & dictionary, const uint & rateCateg) {return dictionary->at(rateCateg)[_dictionaryKey] ;} ;
  uvec * GetInput() {assert(false) ; return NULL;} ;
  void EnterInput(TreeNode * originVertex) {} ;
  void MarkKeyUndefined() {_keyDefined = false ;} ;

  IntermediateNode(): _isSolved(false) {_parent = NULL ;  _keyDefined = false ;};
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   //vec * _solution ;
};

#endif /* INTERMEDIATENODE_H */
