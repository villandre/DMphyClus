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
  virtual void SetSolution() ;
  virtual void InvalidateSolution() ;
  virtual void SetPattern() ;
  virtual void ToggleSolved() ;
  std::vector<longVec> GetSolution() {return _solution ;} ;
  //mapKey GetPattern(uint locusNum) {return _pattern.at(locusNum) ;} ; 
  TreeNode * GetParent() {return _parent ;} ;
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;} ;
  void SetId(uint vertexId) {_id = vertexId ;} ;
  uint GetId() {return _id ;} ;
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
  // std::vector<mapKey> _pattern ;
  std::vector<longVec> _solution ;
};