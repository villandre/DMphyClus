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

typedef std::vector<std::unordered_map<std::size_t, Col<double>>>* solutionDictionaryType ;
typedef std::unordered_map<std::size_t, Col<double>>* nodePatternDictionaryType ;
typedef std::vector<Col<double>> doubleVec ;

class TreeNode
{
public:
  virtual bool IsSolved() = 0;
  virtual bool CanSolve() = 0;
  virtual bool CanFindKey() = 0;
  virtual void AddChild(TreeNode *) = 0 ;
  virtual void RemoveChildren() = 0 ;
  virtual void RemoveChild(TreeNode *) = 0 ;
  virtual void ComputeSolution(solutionDictionaryType &, const mat &, double &, const uint &) = 0 ;
  virtual void InvalidateSolution() = 0;
  virtual void SetSolved(bool) = 0;
  virtual void SetInput(std::vector<uvec> *) = 0 ;
  virtual std::vector<TreeNode *> GetChildren() = 0;
  virtual void DeriveKey(solutionDictionaryType &, const uint &, const uint &, const uint &) = 0;
  virtual vec GetSolution(solutionDictionaryType &, const uint &, const std::size_t &) = 0;
  //virtual void EnterInput(TreeNode *) = 0;
  virtual std::vector<uvec> * GetInput() = 0 ;
  virtual void MarkKeyUndefined() = 0 ;
  
  std::vector<std::size_t> GetDictionaryKeyVec() const { return _dictionaryKeyVec ;}
  TreeNode * GetParent() {return _parent ;}
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;}
  void SetId(uint vertexId) {_id = vertexId ;}
  uint GetId() {return _id ;}
  bool IsKeyDefined() {return _keyDefined ;}
  
  bool GetWithinParentBranch() {return _withinParentBranch ;}
  
  // void EnterCommonInfo(TreeNode * originVertex)
  // {
  //   _withinParentBranch = originVertex->GetWithinParentBranch() ;
  //   _dictionaryKey = originVertex->GetDictionaryKey() ;
  //   _keyDefined = true ;
  // }
  void SetWithinParentBranch(bool parentBranchWithin) {_withinParentBranch = parentBranchWithin ;}
  void DeriveKeys(solutionDictionaryType & solutionDictionary, const uint & matListIndex, const uint & numElements)
  {
    _dictionaryKeyVec.reserve(numElements) ;
    uint rateCategIndex = 0 ;
    for (uint i = 0; i < numElements ; i++)
    {
      DeriveKey(solutionDictionary, rateCategIndex, matListIndex, i) ;
      rateCategIndex = littleCycle(rateCategIndex+1, solutionDictionary->size()) ;
    }
  }
  virtual ~TreeNode() { }
  
  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  bool _withinParentBranch ; // true if the parent branch has within-cluster transition probabilities.
  std::vector<std::size_t> _dictionaryKeyVec ;
  bool _keyDefined ;
};

#endif /* TREENODE_H */
