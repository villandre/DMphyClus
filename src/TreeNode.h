#include <boost/functional/hash.hpp>
#include <assert.h>
#include <RcppArmadillo.h>
#include <unordered_map>

using namespace arma ;

#ifndef TREENODE_H
#define TREENODE_H

typedef std::unordered_map<std::size_t, Col<long double>, std::size_t> solutionDictionaryType ;
typedef std::unordered_map<std::size_t, Col<long double>, std::size_t> nodePatternDictionaryType ;
typedef std::vector<Col<long double>> longVec ;

class TreeNode
{
public:
  virtual bool IsSolved() = 0;
  virtual bool CanSolve() ;
  virtual void AddChild(TreeNode* child) ;
  virtual void RemoveChild(TreeNode* child) ;
  virtual void SetSolution(Col<long double>) ;
  virtual void ComputeSolution() ;
  virtual void InvalidateSolution() ;
  virtual void SetPattern() ;
  virtual void ToggleSolved() ;
  virtual void SetInput(const uvec &) ;
  virtual std::vector<TreeNode *> GetChildren() ;
  virtual void DeriveKey(solutionDictionaryType &) ;
  std::size_t GetDictionaryKey() const { return _dictionaryKey ;};
  virtual Col<long double> GetSolution() ;
  TreeNode * GetParent() {return _parent ;} ;
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;} ;
  void SetId(uint vertexId) {_id = vertexId ;} ;
  uint GetId() {return _id ;} ;
  void SetTransProbMatrix(const mat & transProbMatrix, std::size_t rateCategory, bool withinCluster) {_transProbMatrix = transProbMatrix ; _rateCategory = rateCategory ; _withinCluster = withinCluster ;} ;
  mat GetTransMatrix() {return _transProbMatrix ;} ;

  protected:
    
  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  mat _transProbMatrix ; // This matrix is associated with the supporting branch.
  std::size_t _rateCategory ; // transProbMatrix gives that indication too, but it's easier to have it mentioned explicitly.
  bool _withinCluster ; // Like before, used to hash the pattern, will be converted to std::size_t.
  std::size_t _dictionaryKey ;
};

#endif /* TREENODE_H */