#include <boost/functional/hash.hpp>
#include <RcppArmadillo.h>

using namespace arma ;

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
  virtual std::size_t GetPattern() ;
  virtual Col<long double> GetSolution() ;
  //mapKey GetPattern(uint locusNum) {return _pattern.at(locusNum) ;} ; 
  TreeNode * GetParent() {return _parent ;} ;
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;} ;
  void SetId(uint vertexId) {_id = vertexId ;} ;
  uint GetId() {return _id ;} ;
  void SetTransProbMatrix (const mat & transProbMatrix) { _transProbMatrix = transProbMatrix ;} ;
  mat GetTransMatrix() {return _transProbMatrix ;} ;
  // void InvalidatePattern() 
  // {
  //   _pattern.clear() ;
  //   if (_parent != NULL) {
  //     _parent->InvalidatePattern() ;
  //   }
  // }; 

  protected:
    
  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  mat _transProbMatrix ; // This matrix is associated with the supporting branch.
};