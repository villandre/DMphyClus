#include <gsl/gsl_sf_gamma.h>
#include <boost/functional/hash.hpp>
#include "IntermediateNode.h"

void IntermediateNode::InvalidateSolution()
{
  _isSolved = false ;
 
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

void IntermediateNode::ComputeSolutions(solutionDictionaryType & solutionDictionary, const std::vector<mat> & transProbMats, const uint & transMatIndex)
{
  uint rateCateg = 0 ;
  _exponentIncrementVec.fill(0) ; // This must be reset to 0, since AugTree persists between iterations.
  for (uint i = 0 ; i < _dictionaryIterVec.size() ; i++)
  {
    ComputeSolution(solutionDictionary, transProbMats.at(rateCateg), rateCateg, i, transMatIndex) ;
    rateCateg = littleCycle(rateCateg + 1, transProbMats.size()) ;
  }
  _isSolved = true ;
  _updateFlag = true ;
}

// I'm using a scaling strategy to avoid computational zeros, where when the maximum value in my solution vector
// gets too small, I factorize it out and increment _exponentContainer, whose total value is taken into
// account when computing the likelihood in Forest::ComputeLikelihood.
// Under this strategy, some elements of the L vector may take value 0 before the scaling is applied, 
// but only when they're much smaller than the maximum, in which case, they won't affect the mean significantly.
void IntermediateNode::ComputeSolution(solutionDictionaryType & solutionDictionary, const mat & transProbM, const uint & rateCategory, const uint & elementNum, const uint & transMatrixIndex)
{
  cout << "Entered ComputeSolution! \n" ;
  float exponentIncrement = 0 ;
  uint rateCateg = littleCycle(elementNum, solutionDictionary->size()) ;
  S newS = GetSfromVertex(elementNum, transMatrixIndex, solutionDictionary->size()) ;
  cout << "Looking for solution... " ;
  mapIterator solutionIter = solutionDictionary->at(rateCateg).find(newS) ;
  cout << "Done! \n" ;
  if (solutionIter != solutionDictionary->at(rateCateg).end()) 
  {
    cout << "Computing distance... " ;
    _iterMove.at(elementNum) = std::distance(solutionIter, _dictionaryIterVec.at(elementNum)) ;
    cout << "Done! \n" ;
    if (elementNum < 30) 
    {
      cout << "(Found solution) _iterMove is: " << _iterMove.at(elementNum) << "\n" ;
    }
    _dictionaryIterVec.at(elementNum) = solutionIter ;
    _exponentIncrementVec.at(elementNum) = solutionIter->second.second ;
  }
  else
  {
    vec mySolution(transProbM.n_rows, fill::ones) ;
    
    for (auto & child : _children) 
    {
      mySolution = mySolution % (transProbM*child->GetSolution(elementNum, solutionDictionary->size())) ;
    }
    
    double myMax = max(mySolution) ;
    bool status = myMax < 1e-150 ; // To account for computational zeros... Will only work with bifurcating trees though.
    if (status)
    {
      mySolution = mySolution/myMax ;
      exponentIncrement = log(myMax) ;
      _exponentIncrementVec.at(elementNum) = exponentIncrement;
    }
    std::pair<mapIterator, bool> insertResult = solutionDictionary->at(rateCategory).insert(std::pair<S,std::pair<vec,float>>(GetSfromVertex(elementNum, transMatrixIndex, solutionDictionary->size()), std::pair<vec, float>(mySolution, exponentIncrement))) ;
    _iterMove.at(elementNum) = std::distance(insertResult.first, _dictionaryIterVec.at(elementNum)) ;
    if (elementNum < 30) 
    {
      cout << "_iterMove is: " << _iterMove.at(elementNum) << "\n" ;
    }
    _dictionaryIterVec.at(elementNum) = insertResult.first ;
    
  }
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

S IntermediateNode::GetSfromVertex(const uint & elementNum, const uint & transMatIndex, const uint & numRateCats)
{
  bool childrenWithinCluster = _children.at(0)->GetWithinParentBranch() ;
  
  return S(std::hash<S>{} (_children.at(0)->GetDictionaryIterator(elementNum, numRateCats)->first),
           std::hash<S>{} (_children.at(1)->GetDictionaryIterator(elementNum, numRateCats)->first),
           childrenWithinCluster,
           transMatIndex) ;
}
