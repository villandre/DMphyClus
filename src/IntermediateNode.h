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
  void ComputeSolution(solutionDictionaryType &, const mat &, double &, const uint &, const uint &) ;
  void InvalidateSolution() ;
  void SetSolved(bool status) {_isSolved = status ;};
  void SetInput(std::vector<uvec> *) { assert(false) ;};
  std::vector<TreeNode *> GetChildren() {return _children;};
  void DeriveKey(solutionDictionaryType &, const uint &, const uint &, const uint &) ;
  vec GetSolution(solutionDictionaryType & dictionary, const uint & rateCateg, const std::size_t & dictionaryKey) {return dictionary->at(rateCateg)[dictionaryKey] ;} ;
  std::vector<uvec> * GetInput() {assert(false) ; return NULL;} ;
  void EnterInput(TreeNode * originVertex) {} ;
  void MarkKeyUndefined() {_keyDefined = false ;} ;

  IntermediateNode(uint & numLoci, uint & numRates): _isSolved(false) {
    _parent = NULL ;  
    _keyDefined = false ;
    _dictionaryKeyVec.reserve(numLoci*numRates) ;
  };
  
protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
};

#endif /* INTERMEDIATENODE_H */
