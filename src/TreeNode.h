// [[Rcpp::depends(RcppArmadillo)]]

#include <boost/functional/hash.hpp>
#include <assert.h>
#include <RcppArmadillo.h>
#include <unordered_map>
#include <gsl/gsl_rng.h>

using namespace arma ;

#ifndef TREENODE_H
#define TREENODE_H

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
  //virtual void SetSolution(vec &) = 0 ;
  virtual void ComputeSolution(solutionDictionaryType &, const mat &, double *) = 0 ;
  virtual void InvalidateSolution() = 0;
  virtual void SetSolved(bool) = 0;
  virtual void SetInput(uvec *) = 0 ;
  virtual std::vector<TreeNode *> GetChildren() = 0;
  virtual void DeriveKey(solutionDictionaryType &) = 0;
  virtual vec GetSolution(solutionDictionaryType &) = 0;
  virtual void EnterInput(TreeNode *) = 0;
  virtual uvec * GetInput() = 0 ;
  
  std::size_t GetDictionaryKey() const { return _dictionaryKey ;};
  TreeNode * GetParent() {return _parent ;} ;
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;} ;
  void SetId(uint vertexId) {_id = vertexId ;} ;
  uint GetId() {return _id ;} ;
  bool IsKeyDefined() {return _keyDefined ;} ;
  bool GetWithinParentBranch() {return _withinParentBranch ;} ;
  std::size_t GetRateCategory() {return _rateCategory ;} ;
  void MarkKeyUndefined() { _keyDefined = false ;} ;
  
  void EnterCommonInfo(TreeNode * originVertex)
  {
    _id = originVertex->GetId() ;
    _rateCategory = originVertex->GetRateCategory() ;
    _withinParentBranch = originVertex->GetWithinParentBranch() ;
    _dictionaryKey = originVertex->GetDictionaryKey() ;
    _keyDefined = true ;
  }
  void SetWithinParentBranchAndRateCateg(bool parentBranchWithin, uint rateCategory) {_withinParentBranch = parentBranchWithin ; _rateCategory = rateCategory ;} ;
  void SetWithinParentBranch(bool parentBranchWithin) {_withinParentBranch = parentBranchWithin ;} ;
  virtual ~TreeNode() { };
  
  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  std::size_t _rateCategory ; // transProbMatrix gives that indication too, but it's easier to have it mentioned explicitly. This info is used to hash nodes.
  bool _withinParentBranch ; // true if the parent branch has within-cluster transition probabilities.
  std::size_t _dictionaryKey ;
  bool _keyDefined ;
};

#endif /* TREENODE_H */
