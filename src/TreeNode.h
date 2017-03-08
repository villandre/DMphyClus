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

typedef std::vector<std::map<std::size_t, Col<double>>>* solutionDictionaryType ;
typedef std::vector<Col<double>> doubleVec ;
typedef std::vector<std::map<std::size_t, vec>::iterator> iterVec ;
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
  //virtual void DeriveKey(solutionDictionaryType &, const uint &, const uint &, const uint &) = 0;
  virtual vec GetSolution(solutionDictionaryType &, const uint &, const std::size_t &) = 0;
  //virtual void EnterInput(TreeNode *) = 0;
  virtual std::vector<uvec> * GetInput() = 0 ;
  virtual void MarkKeyUndefined() = 0 ;
  virtual void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, vec &) = 0 ;
  virtual void UpdateMapAndIterVec(solutionDictionaryType &, uint &) = 0;
  
  iterVec GetDictionaryIterVec() { return _dictionaryIterVec ;}
  TreeNode * GetParent() {return _parent ;}
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;}
  void SetId(uint vertexId) {_id = vertexId ;}
  uint GetId() {return _id ;}
  bool IsKeyDefined() {return _keyDefined ;}
  
  bool GetWithinParentBranch() {return _withinParentBranch ;}
  
  void SetWithinParentBranch(bool parentBranchWithin) {_withinParentBranch = parentBranchWithin ;}
  // void DeriveKeys(solutionDictionaryType & solutionDictionary, const uint & matListIndex, const uint & numElements)
  // {
  //   _dictionaryIterVec.reserve(numElements) ;
  //   _iteratorMove.reserve(numElements) ;
  //   uint rateCategIndex = 0 ;
  //   for (uint i = 0; i < numElements ; i++)
  //   {
  //     DeriveKey(solutionDictionary, rateCategIndex, matListIndex, i) ;
  //     rateCategIndex = littleCycle(rateCategIndex+1, solutionDictionary->size()) ;
  //   }
  // }
  virtual ~TreeNode() { }
  
  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  bool _withinParentBranch ; // true if the parent branch has within-cluster transition probabilities.
  iterVec _dictionaryIterVec ;
  //bool _keyDefined ;
  std::vector<unsigned int> _iteratorMove ;
};

#endif /* TREENODE_H */
