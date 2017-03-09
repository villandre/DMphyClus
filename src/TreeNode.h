// [[Rcpp::depends(RcppArmadillo)]]

#include <boost/functional/hash.hpp>
#include <assert.h>
#include <RcppArmadillo.h>
#include <unordered_map>
#include <gsl/gsl_rng.h>

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
};

bool operator==(const S& lhs, const S& rhs) 
{
  return (((lhs._aHash == rhs._aHash) && (lhs._otherHash == rhs._otherHash)) || 
          ((lhs._aHash == rhs._otherHash) && (lhs._otherHash == rhs._aHash))) &&
          (lhs._withinFlag == rhs._withinFlag) && (lhs._transMatrixIndex == rhs._transMatrixIndex) ;
}

struct classcomp {
  bool operator() (const S& lhs, const S& rhs) const
  {
    if (lhs._transMatrixIndex != rhs._transMatrixIndex) return true ;
    if (lhs._withinFlag != rhs._withinFlag) return true ;
    if (!(((lhs._aHash == rhs._aHash)&&(lhs._otherHash == rhs._otherHash))||((lhs._aHash == rhs._otherHash)&&(lhs._otherHash == rhs._aHash)))) return true ;
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
    std::vector<std::size_t> hashKeys ;
    hashKeys.reserve(4) ;
    std::vector<std::size_t> allInfoForHashing{s._aHash,s._otherHash,s._transMatrixIndex,s._withinFlag} ;
    
    std::sort(allInfoForHashing.begin(), allInfoForHashing.begin()+2); // The children keys should be re-ordered, not the within-cluster indicator. This will make the function symmetrical.
    return boost::hash_range(allInfoForHashing.begin(), allInfoForHashing.end()) ;
  };
} ;
}

typedef std::vector<std::unordered_map<S, vec>>* solutionDictionaryType ;
typedef std::vector<vec> doubleVec ;
typedef std::unordered_map<S, vec>::iterator mapIterator ;
typedef std::vector<mapIterator> iterVec ;

class TreeNode
{
public:
  virtual bool IsSolved() = 0;
  virtual bool CanSolve() = 0;
  virtual void AddChild(TreeNode *) = 0 ;
  virtual void RemoveChildren() = 0 ;
  virtual void RemoveChild(TreeNode *) = 0 ;
  virtual void ComputeSolution(solutionDictionaryType &, const mat &, double &, const uint &, const uint &, const uint &) = 0 ;
  virtual void InvalidateSolution() = 0;
  virtual void SetSolved(bool) = 0;
  virtual void SetInput(std::vector<uvec> *) = 0 ;
  virtual std::vector<TreeNode *> GetChildren() = 0;
  virtual vec GetSolution(const uint &, const uint & numRateCats) = 0;
  
  virtual std::vector<uvec> * GetInput() = 0 ;
  virtual void MarkKeyUndefined() = 0 ;
  virtual void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, vec &, const uint &) = 0 ;
  virtual void InitMapAndIterVec(solutionDictionaryType &) = 0;
  virtual void CopyIterVec() = 0 ;
  virtual void RestoreIterVec() = 0;
  
  std::vector<std::pair<const S, vec>*> GetDictionaryIterVec() { return _dictionaryIterVec ;}
  TreeNode * GetParent() {return _parent ;}
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;}
  void SetId(uint vertexId) {_id = vertexId ;}
  uint GetId() {return _id ;}
  void NegateFlag() {_updateFlag = false ;} 
  
  bool GetWithinParentBranch() {return _withinParentBranch ;}
  
  void SetWithinParentBranch(bool parentBranchWithin) {_withinParentBranch = parentBranchWithin ;}
  
  virtual ~TreeNode() { }
  
  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  bool _withinParentBranch ; // true if the parent branch has within-cluster transition probabilities.
  std::vector<std::pair<const S, vec>*> _dictionaryIterVec ;
  std::vector<std::pair<const S, vec>*> _previousIterVec ;
  bool _updateFlag ;
};

#endif /* TREENODE_H */
