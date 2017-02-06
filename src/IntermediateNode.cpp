#include "IntermediateNode.h"
#include <RcppArmadillo.h>

void IntermediateNode::InvalidateSolution() 
{
  _isSolved = FALSE ;
  _solution = zeros<Col<long double>>(_solution.size()) ;
  if (_parent != NULL) { // Root has a NULL parent. 
    _parent->InvalidateSolution() ; 
  }
}

bool IntermediateNode::CanSolve() {
  
  std::vector<bool> childDefined(_children.size()) ; 
  std::transform(_children.begin(), _children.end(), childDefined.begin(), [] (TreeNode * childNodePointer) {return childNodePointer->IsSolved() ;}) ; 
  return std::all_of(childDefined.begin(), childDefined.end(), [](bool v) { return v; });
}

void IntermediateNode::AddChild(TreeNode* child) {
  
  _children.push_back(child) ;
}

void IntermediateNode::RemoveChild(TreeNode* child) 
{
  auto childPos = find(_children.begin(), _children.end(), child);
  _children.erase(childPos) ;
}

void IntermediateNode::SetSolution(mat & transProbMatrix, std::unordered_map<size_t, Col<long double>> & dictionary) 
{
  if (dictionary.find(_pattern) != dictionary.end()) {
    _solution = dictionary[_pattern] ;
  } 
  std::vector<Col<long double>> childrenCombined(_children.size()) ;
  std::transform(_children.begin(), _children.end(), childrenCombined, [&transProbMatrix] (TreeNode * child) {return transProbMatrix*child->GetSolution() ;}) ;  
  _solution = accumulate(childrenCombined.begin(), childrenCombined.end(), childrenCombined.at(0), [] (Col<long double> & a, Col<long double> & b) {return a%b ;}) ;
}

void IntermediateNode::SetPattern(std::unordered_map<size_t, Col<long double>> & dictionary) 
{
  std::vector<int> childrenHashes(_children.size()) ;
  //TO_DO
}