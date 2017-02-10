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

void IntermediateNode::ComputeSolution() 
{
  Col<long double> mySolution(_transProbMatrix.n_rows, fill::ones) ;
  for(auto & child : _children) 
  {
    mySolution = mySolution % child->GetTransMatrix()*child->GetSolution() ;
  }
  _solution = mySolution ;
}

void IntermediateNode::SetPattern(std::unordered_map<mapKey, Col<long double>, MyHash> & dictionary) 
{
  std::vector<int> childrenHashes(_children.size()) ;
  //TO_DO
}