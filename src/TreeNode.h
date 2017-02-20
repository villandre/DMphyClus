// [[Rcpp::depends(RcppArmadillo)]]

#include <boost/functional/hash.hpp>
#include <assert.h>
#include <RcppArmadillo.h>
#include <unordered_map>

using namespace arma ;

#ifndef TREENODE_H
#define TREENODE_H

typedef std::unordered_map<std::size_t, Col<double>> solutionDictionaryType ;
typedef std::unordered_map<std::size_t, Col<double>> nodePatternDictionaryType ;
typedef std::vector<Col<double>> doubleVec ;

class TreeNode
{
public:
  virtual bool IsSolved() = 0;
  virtual bool CanSolve() = 0;
  virtual bool CanFindKey() = 0;
  virtual void AddChild(TreeNode *) = 0 ;
  virtual void RemoveChildren() = 0 ;
  virtual void SetSolution(vec &) = 0 ;
  virtual void ComputeSolution(solutionDictionaryType &) = 0 ;
  virtual void InvalidateSolution() = 0;
  virtual void ToggleSolved() = 0;
  virtual void SetInput(const uvec &) = 0 ;
  virtual std::vector<TreeNode *> GetChildren() = 0;
  virtual void DeriveKey(solutionDictionaryType &) = 0;
  virtual vec GetSolution() = 0;

  std::size_t GetDictionaryKey() const { return _dictionaryKey ;};
  TreeNode * GetParent() {return _parent ;} ;
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;} ;
  void SetId(uint vertexId) {_id = vertexId ;} ;
  uint GetId() {return _id ;} ;
  void SetTransProbMatrix(const mat & transProbMatrix, std::size_t rateCategory, bool withinCluster) {_transProbMatrix = transProbMatrix ; _rateCategory = rateCategory ; _withinParentBranch = withinCluster ;} ;
  mat GetTransMatrix() {return _transProbMatrix ;} ;
  bool IsKeyDefined() {return _keyDefined ;} ;
  bool GetWithinParentBranch() {return _withinParentBranch ;} ;
  
  virtual ~TreeNode() { };
  
  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  mat _transProbMatrix ; // This matrix is associated with the supporting branch.
  std::size_t _rateCategory ; // transProbMatrix gives that indication too, but it's easier to have it mentioned explicitly.
  bool _withinParentBranch ; // true if the parent branch has within-cluster transition probabilities.
  std::size_t _dictionaryKey ;
  bool _keyDefined ;
};

#endif /* TREENODE_H */
