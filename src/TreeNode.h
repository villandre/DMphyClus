#include <boost/functional/hash.hpp>
#include <RcppArmadillo.h>

using namespace arma ;

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
  Col<long double> GetSolution() {return _solution ;} ;
  std::size_t GetPattern() {return _pattern;} ; 
  TreeNode * GetParent() {return _parent ;} ;
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;} ;
  void SetId(uint vertexId) {_id = vertexId ;} ;
  void AssignSolution(const Col<long double> & knownSolution) {_solution = knownSolution ;} ;
  void InvalidatePattern() 
  {
    _pattern = 0 ;
    _parent->InvalidatePattern() ;
  }; 
  
  std::vector<TreeNode *> _children ;

  protected:
    
  uint _id ; // From 1 to number of nodes. Used for exporting the phylogeny to R.
  TreeNode * _parent ;
  size_t _pattern ;
  Col<long double> _solution ;
};