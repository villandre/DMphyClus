#include <gsl/gsl_sf_gamma.h>
#include <boost/functional/hash.hpp>
#include "IntermediateNode.h"

void IntermediateNode::InvalidateSolution()
{
  _isSolved = false ;
  _solution = zeros<Col<double>>(_solution.size()) ;
  if (_parent != NULL)
  { // Root has a NULL parent.
    _parent->InvalidateSolution() ;
  }
}

bool IntermediateNode::CanSolve()
{
  std::vector<bool> childDefined(_children.size()) ;
  childDefined.reserve(_children.size()) ;
  for (auto & i : _children)
  {
    childDefined.push_back(i->IsSolved()) ;
  }
  return std::all_of(childDefined.begin(), childDefined.end(), [](bool v) { return v; });
}

bool IntermediateNode::CanFindKey()
{
  std::vector<bool> childKeyDefined(_children.size()) ;
  childKeyDefined.reserve(_children.size()) ;
  for (auto & i : _children)
  {
    childKeyDefined.push_back(i->IsKeyDefined()) ; 
  }
  return std::all_of(childKeyDefined.begin(), childKeyDefined.end(), [](bool v) { return v; });
}

// I'm using a scaling strategy to avoid computational zeros, where when the maximum value in my solution vector
// gets too small, I factorize it out and increment _exponentContainer, whose total value is taken into
// account when computing the likelihood in Forest::ComputeLikelihood.
// Under this strategy, some elements of the L vector may take value 0, but only when they're much smaller
// than the maximum, in which case, they won't affect the mean significantly.
void IntermediateNode::ComputeSolution(solutionDictionaryType & solutionDictionary, const mat & transProbM, double * expContainer)
{
  Col<double> mySolution(transProbM.n_rows, fill::ones) ;
  for(auto & child : _children)
  {
    mySolution = mySolution % (transProbM*child->GetSolution()) ;
  }
  double myMax = max(mySolution) ;
  bool status = myMax < 1e-10 ; // To account for computational zeros...
  if (status)
  {
    mySolution = mySolution/myMax ;
    *expContainer = *expContainer + log(myMax);
  }
  _solution = mySolution ;
  _isSolved = true ;
  //solutionDictionary[_dictionaryKey] = mySolution ;
}

void IntermediateNode::DeriveKey(solutionDictionaryType & solutionDictionary)
{
  std::vector<std::size_t> permutations(gsl_sf_fact(_children.size())+2) ; // The hash key is computed from the rate category and the within/between status, hence +2.
  std::vector<std::size_t> hashKeys ;
  hashKeys.reserve(permutations.size()) ;
  std::transform(_children.begin(), _children.end(), permutations.begin(), [] (TreeNode * childPointer) {return childPointer->GetDictionaryKey() ;}) ;
  permutations[_children.size()] = _rateCategory ;
  permutations[_children.size() + 1] = (std::size_t) _withinParentBranch ;
  std::sort(permutations.begin(), permutations.end()-2); // The children keys should be re-ordered, not the rate category index and within-cluster indicator.
  do {
    hashKeys.push_back(boost::hash_range(permutations.begin(), permutations.end())) ;
  } while ( std::next_permutation(permutations.begin(),permutations.end() - 2) );
  bool foundSolution = false ;
  for(auto & hashKey : hashKeys) {
    if (solutionDictionary->find(hashKey) != solutionDictionary->end()) {
      _dictionaryKey = hashKey ;
      foundSolution = true ;
      break ;
    }
  }
  if (!foundSolution) {
    _dictionaryKey = hashKeys[0] ;
  }
  _keyDefined = true ;
}

void IntermediateNode::RemoveChild(TreeNode * childToRemove)
{
  auto ChildIterPos = std::find(_children.begin(), _children.end(), childToRemove) ;
  if (ChildIterPos == _children.end())
  {
    cerr << "Warning: Trying to remove a child that was not found! \n" ;
  }
  else 
  {
    _children.erase(ChildIterPos) ;
  }
}