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
  std::copy(_dictionaryIterVec.begin(), _dictionaryIterVec.end(), _previousIterVec.begin()) ;
  for (uint i = 0 ; i < _dictionaryIterVec.size() ; i++)
  {
    ComputeSolution(solutionDictionary, transProbMats, i, transMatIndex) ;
  }
  _isSolved = true ;
  _updateFlag = true ;
}

// I'm using a scaling strategy to avoid computational zeros, where when the maximum value in my solution vector
// gets too small, I factorize it out and increment _exponentContainer, whose total value is taken into
// account when computing the likelihood in Forest::ComputeLikelihood.
// Under this strategy, some elements of the L vector may take value 0 before the scaling is applied, 
// but only when they're much smaller than the maximum, in which case, they won't affect the mean significantly.
void IntermediateNode::ComputeSolution(solutionDictionaryType & solutionDictionary, const std::vector<mat> & transProbMatVec, const uint & locusNum, const uint & transMatrixIndex)
{
  S newS = GetSfromVertex(locusNum, transMatrixIndex, transProbMatVec.size()) ;
  mapIterator solutionIter = solutionDictionary->find(newS) ;
  if (solutionIter != solutionDictionary->end()) 
  {
    _dictionaryIterVec.at(locusNum) = solutionIter ;
    // if ((_id == 6) && (locusNum < 10))
    // {
    //   solutionIter->second.at(0).first.print("Found solution:") ;
    // }
  }
  else
  {
    std::vector<std::pair<vec, float>> mySolution(transProbMatVec.size(), std::pair<vec,float>(vec(transProbMatVec.at(0).n_rows, fill::ones),0)) ;
    for (auto & child : _children) 
    {
      for (uint rateIndex = 0 ; rateIndex < mySolution.size() ; rateIndex++)
      {
        mySolution.at(rateIndex).first = mySolution.at(rateIndex).first % (transProbMatVec.at(rateIndex)*child->GetSolution(locusNum, rateIndex)) ;
        double myMax = max(mySolution.at(rateIndex).first) ;
        bool status = myMax < 1e-150 ; // To account for computational zeros... Will only work with bifurcating trees though.
        if (status)
        {
          mySolution.at(rateIndex).first = mySolution.at(rateIndex).first/myMax ;
          mySolution.at(rateIndex).second = log(myMax) ;
        }
      }
    }
    std::pair<mapIterator, bool> insertResult = solutionDictionary->insert(std::pair<S,mapContentType>(GetSfromVertex(locusNum, transMatrixIndex, transProbMatVec.size()), mapContentType(mySolution))) ;
    
    _dictionaryIterVec.at(locusNum) = insertResult.first ;
    // if ((_id == 6) && (locusNum < 10))
    // {
    // _dictionaryIterVec.at(locusNum)->second.at(0).first.print("New solution:") ;
    // }
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
