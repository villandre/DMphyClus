#include <gsl/gsl_sf_gamma.h>
#include <boost/functional/hash.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/future.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread/lock_guard.hpp> 
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

void IntermediateNode::ComputeSolutions(solutionDictionaryType & solutionDictionary, const std::vector<mat> & transProbMats, const uint & transMatIndex, boost::asio::io_service & ioService, boost::mutex & myMutex)
{
  cout << "Computing solutions for node " << _id << endl ; 
  std::copy(_dictionaryIterVec.begin(), _dictionaryIterVec.end(), _previousIterVec.begin()) ;
  std::vector<boost::unique_future<bool>> pending_data;
  for (uint i = 0 ; i < _dictionaryIterVec.size() ; i++)
  {
    typedef boost::packaged_task<bool> task_t;
    
    boost::shared_ptr<task_t> task = boost::make_shared<task_t>(
      boost::bind(&IntermediateNode::ComputeSolution, this, solutionDictionary, std::cref(transProbMats), std::cref(i), std::cref(transMatIndex), std::ref(myMutex)));
    
    boost::unique_future<bool> fut = task->get_future();
    
    pending_data.push_back(std::move(fut));
    ioService.post(boost::bind(&task_t::operator(), task));
    
    //ioService.post(boost::bind(&IntermediateNode::ComputeSolution, this, solutionDictionary, std::cref(transProbMats), std::cref(i), std::cref(transMatIndex), std::ref(myMutex)));
    //ComputeSolution(solutionDictionary, transProbMats, i, transMatIndex, myMutex) ;
  }
  boost::wait_for_all(pending_data.begin(), pending_data.end()); 
  cout << "Done! \n" ;
  _isSolved = true ;
  _updateFlag = true ;
}

// I'm using a scaling strategy to avoid computational zeros, where when the maximum value in my solution vector
// gets too small, I factorize it out and increment _exponentContainer, whose total value is taken into
// account when computing the likelihood in Forest::ComputeLikelihood.
// Under this strategy, some elements of the L vector may take value 0 before the scaling is applied, 
// but only when they're much smaller than the maximum, in which case, they won't affect the mean significantly.
bool IntermediateNode::ComputeSolution(solutionDictionaryType & solutionDictionary, const std::vector<mat> & transProbMatVec, const uint & locusNum, const uint & transMatrixIndex, boost::mutex & myMutex)
{
  S newS = GetSfromVertex(locusNum, transMatrixIndex, transProbMatVec.size()) ;
  mapIterator solutionIter ;
  {
    boost::lock_guard<boost::mutex> lock(myMutex);
    solutionIter = solutionDictionary->find(newS) ;
  }
  if (solutionIter != solutionDictionary->end()) 
  {
    boost::lock_guard<boost::mutex> lock(myMutex);
    _dictionaryIterVec.at(locusNum) = solutionIter ;
  }
  else
  {
    std::vector<std::pair<vec, float>> mySolution(transProbMatVec.size(), std::pair<vec,float>(vec(transProbMatVec.at(0).n_rows, fill::ones),0)) ;
    for (auto & child : _children) 
    {
      for (uint rateIndex = 0 ; rateIndex < mySolution.size() ; rateIndex++)
      {
        {
          boost::lock_guard<boost::mutex> lock(myMutex);
          mySolution.at(rateIndex).first = mySolution.at(rateIndex).first % (transProbMatVec.at(rateIndex)*child->GetSolution(locusNum, rateIndex, myMutex)) ;
        }
        double myMax = max(mySolution.at(rateIndex).first) ;
        bool status = myMax < 1e-150 ; // To account for computational zeros... Will only work with bifurcating trees though.
        if (status)
        {
          mySolution.at(rateIndex).first = mySolution.at(rateIndex).first/myMax ;
          mySolution.at(rateIndex).second = log(myMax) ;
        }
      }
    }
    {
      boost::lock_guard<boost::mutex> lock(myMutex);
      std::pair<mapIterator, bool> insertResult = solutionDictionary->insert(std::pair<S,mapContentType>(GetSfromVertex(locusNum, transMatrixIndex, transProbMatVec.size()), mapContentType(mySolution))) ;
      _dictionaryIterVec.at(locusNum) = insertResult.first ;
    }
  }
  return true ;
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
