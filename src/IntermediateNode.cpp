#include <RcppArmadillo.h>
#include <gsl/gsl_sf_gamma.h>
#include "IntermediateNode.h"

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

void IntermediateNode::DeriveKey(solutionDictionaryType & solutionDictionary)
{
  std::vector<std::size_t> permutations(gsl_sf_fact(_children.size())) ;
  std::vector<std::size_t> hashKeys ;
  hashKeys.reserve(permutations.size()) ;
  std::transform(_children.begin(), _children.end(), permutations.begin(), [] (TreeNode * childPointer) {return childPointer->GetDictionaryKey() ;}) ;
  permutations.at(_children.size()) = _rateCategory ;
  permutations.at(_children.size() + 1) = (std::size_t) _withinCluster ;
  std::sort(permutations.begin(), permutations.end()-2); // The children keys should be re-ordered, not the rate category index and within-cluster indicator.
  do {
    hashKeys.push_back(std::hash<std::vector<size_t>>(permutations)) ;
  } while ( std::next_permutation(permutations.begin(),permutations.end() - 2) );
  bool foundSolution = false ;
  for(auto & hashKey : hashKeys) {
    if (solutionDictionary.find(hashKey) != solutionDictionary.end()) {
      _dictionaryKey = hashKey ;
      foundSolution = true ;
      break ;
    }
  }
  if (!foundSolution) {
    _dictionaryKey = hashKeys.at(0) ;
  }
}