// [[Rcpp::depends(RcppArmadillo)]]

#include <boost/functional/hash.hpp>
#include <assert.h>
#include <RcppArmadillo.h>
#include <unordered_map>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <boost/bind.hpp>
#include "threadpool.h"

using namespace arma ;

#ifndef TREENODE_H
#define TREENODE_H

#include "helper.h"

struct S {
  std::size_t _aHash;
  std::size_t _otherHash;
  std::size_t _withinFlag;
  std::size_t _transMatrixIndex;
  
  S(std::size_t aHash, std::size_t otherHash, bool withinFlag, uint transMatrixIndex) : _aHash(aHash), _otherHash(otherHash) 
  {
    _withinFlag = (std::size_t) withinFlag ;
    _transMatrixIndex = (std::size_t) transMatrixIndex ;
  }
  S() {};
  void print() const
  {
    cout << "S elements: " << _aHash << ", " << _otherHash << ", " << _transMatrixIndex << ", " << _withinFlag << ".\n" ;
  }
};

struct classcomp {
  bool operator() (const S& lhs, const S& rhs) const
  {
    if (lhs._transMatrixIndex != rhs._transMatrixIndex) return lhs._transMatrixIndex < rhs._transMatrixIndex ;
    if (lhs._withinFlag != rhs._withinFlag) return lhs._withinFlag < rhs._withinFlag ;
    if (std::min(lhs._aHash, lhs._otherHash) != std::min(rhs._aHash, rhs._otherHash)) return std::min(lhs._aHash,lhs._otherHash) < std::min(rhs._aHash, rhs._otherHash) ;
    if (std::max(lhs._aHash,lhs._otherHash) != std::max(rhs._aHash, rhs._otherHash)) return std::max(lhs._aHash,lhs._otherHash) < std::max(rhs._aHash, rhs._otherHash) ;
    return false ;
  }
};

// custom specialization of std::hash can be injected in namespace std
namespace std
{
template<> struct hash<S>
{
  typedef S argument_type;
  typedef std::size_t result_type;
  result_type operator()(argument_type const& s) const
  {
    std::vector<std::size_t> allInfoForHashing{s._aHash,s._otherHash,s._transMatrixIndex,s._withinFlag} ;
    
    std::sort(allInfoForHashing.begin(), allInfoForHashing.begin()+2); // The children keys should be re-ordered, not the within-cluster indicator. This will make the function symmetrical.
    std::size_t seed = 0 ;
    boost::hash_combine(seed, allInfoForHashing.at(0));
    boost::hash_combine(seed, allInfoForHashing.at(1));
    boost::hash_combine(seed, allInfoForHashing.at(2));
    boost::hash_combine(seed, allInfoForHashing.at(3));
    
    return seed;
    //return boost::hash_range(allInfoForHashing.begin(), allInfoForHashing.end()) ;
  };
} ;
}

typedef std::vector<std::pair<vec, float>> mapContentType ;
typedef std::map<S, mapContentType, classcomp>* solutionDictionaryType ;
typedef std::map<S, mapContentType, classcomp> solutionDictionaryTypeNoPoint ;
typedef std::vector<vec> doubleVec ;
typedef std::map<S, mapContentType, classcomp>::const_iterator mapIterator ;
typedef std::vector<mapIterator> iterVec ;


class TreeNode
{
public:
  virtual void AddChild(TreeNode *) = 0 ;
  virtual void RemoveChild(TreeNode *) = 0 ;
  virtual std::vector<TreeNode *> GetChildren() = 0;
  virtual void RemoveChildren() = 0 ;
  
  virtual bool IsSolved() = 0;
  virtual bool CanSolve() = 0;
  virtual void SetSolved(bool) = 0;
  virtual void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, const uint &, threadpool_t *, pthread_spinlock_t &) = 0 ;
  virtual void ComputeSolution(solutionDictionaryType &, const std::vector<mat> &, const uint &, const uint &, pthread_spinlock_t &) = 0 ;
  virtual void InvalidateSolution() = 0;
  
  virtual void SetInput(std::vector<uvec> *) = 0 ;
  virtual std::vector<uvec> * GetInput() = 0 ;
  
  virtual void InitMapAndIterVec(solutionDictionaryType &, const uint &) = 0;
  virtual void RestoreIterVecAndExp() = 0;
  virtual mapIterator GetDictionaryIterator(const uint &, const uint &) = 0 ;
  virtual S GetSfromVertex(const uint &, const uint &, const uint &) = 0;
  
  void ComputeSolutionWrap(void* arguments) 
  {
    auto f1 = std::bind(TreeNode::ComputeSolution, this, );
    ComputeSolution(myArgs->_solutionDictionary, myArgs->_transProbMatVec, myArgs->_locusNum, myArgs->_transMatrixIndex, myArgs->_threadpool) ;
  }
  TreeNode * GetParent() {return _parent ;}
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;}
  void SetId(uint vertexId) {_id = vertexId ;}
  uint GetId() {return _id ;}
  void NegateFlag() {_updateFlag = false ;} 
  vec GetSolution(const uint & locusNum, const uint & rateCat) 
  {
    return _dictionaryIterVec.at(locusNum)->second.at(rateCat).first ;
  } 
  vec GetSolutionNoMutex(const uint & locusNum, const uint & rateCat) 
  {
    return _dictionaryIterVec.at(locusNum)->second.at(rateCat).first ;
  } 
  float GetExponent(const uint & locusNum, const uint & rateCat) { return _dictionaryIterVec.at(locusNum)->second.at(rateCat).second ;}
  
  bool GetWithinParentBranch() {return _withinParentBranch ;}
  
  void SetWithinParentBranch(bool parentBranchWithin) {_withinParentBranch = parentBranchWithin ;}
  
  virtual ~TreeNode() { }
  
  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  bool _withinParentBranch ; // true if the parent branch has within-cluster transition probabilities.
  iterVec _dictionaryIterVec ;
  
  //std::vector<float> _iterMove ;
  bool _updateFlag ;
};

#endif /* TREENODE_H */
