#include <gsl/gsl_sf_gamma.h>
#include <boost/functional/hash.hpp>
#include "IntermediateNode.h"

void IntermediateNode::InvalidateSolution()
{
  _isSolved = false ;
  _keyDefined = false ; // If a solution is invalidated, it means something changed in lower vertice, which also invalidates the key.
  
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
// Under this strategy, some elements of the L vector may take value 0 before the scaling is applied, 
// but only when they're much smallerthan the maximum, in which case, they won't affect the mean significantly.
void IntermediateNode::ComputeSolution(solutionDictionaryType & solutionDictionary, const mat & transProbM, double * expContainer, const uint & rateCategory)
{
  Col<double> mySolution(transProbM.n_rows, fill::ones) ;
  for(auto & child : _children)
  {
    mySolution = mySolution % (transProbM*child->GetSolution(solutionDictionary, rateCategory)) ;
  }
  double myMax = max(mySolution) ;
  bool status = myMax < 1e-150 ; // To account for computational zeros... Will only work with bifurcating trees though.
  if (status)
  {
    mySolution = mySolution/myMax ;
    *expContainer = *expContainer + log(myMax);
  }
  
  (*solutionDictionary).at(rateCategory)[_dictionaryKey] = mySolution ;
  //_solution = &(*solutionDictionary).at(_rateCategory)[_dictionaryKey] ;
  _isSolved = true ;
}

void IntermediateNode::DeriveKey(solutionDictionaryType & solutionDictionary, const uint & rateCategory, const uint & matListIndex)
{
  std::vector<std::size_t> childrenKeysAndWithinBetweenFlag(4) ; // The hash key is also computed from the within/between status, and betweenMatListIndex or withinMatListIndex, hence +2. The tree is assumed bifurcating, hence the 4.
  std::vector<std::size_t> hashKeys ;
  hashKeys.reserve(childrenKeysAndWithinBetweenFlag.size()) ;
  std::transform(_children.begin(), _children.end(), childrenKeysAndWithinBetweenFlag.begin(), [] (TreeNode * childPointer) {return childPointer->GetDictionaryKey() ;}) ;
  childrenKeysAndWithinBetweenFlag.at(2) = (std::size_t) _children.at(0)->GetWithinParentBranch() ;
  childrenKeysAndWithinBetweenFlag.at(3) = matListIndex ;
  
  std::sort(childrenKeysAndWithinBetweenFlag.begin(), childrenKeysAndWithinBetweenFlag.end()-2); // The children keys should be re-ordered, not the within-cluster indicator.
  do {
    hashKeys.push_back(boost::hash_range(childrenKeysAndWithinBetweenFlag.begin(), childrenKeysAndWithinBetweenFlag.end())) ;
  } while ( std::next_permutation(childrenKeysAndWithinBetweenFlag.begin(), childrenKeysAndWithinBetweenFlag.end()-2));
  bool foundSolution = false ;
  for(auto & hashKey : hashKeys) {
    if (solutionDictionary->at(rateCategory).find(hashKey) != solutionDictionary->at(rateCategory).end()) {
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