#include <boost/functional/hash.hpp>
#include <assert.h>
#include "TreeNode.h"

class IntermediateNode:public TreeNode
{
public:
  virtual bool IsSolved() {return _isSolved ;};
  virtual bool CanSolve() ;
  virtual std::size_t GetPattern() ; // size_t * are pointers to the hash keys attributed to the children nodes by boost::hash. Beware of collisions!
  virtual void AddChild(TreeNode*) ;
  virtual void RemoveChild(TreeNode*) ;
  virtual void SetSolution(mat &, std::unordered_map<size_t, Col<long double>> &) ;
  virtual void InvalidateSolution() ;
  virtual void SetPattern(std::unordered_map<size_t, Col<long double>> &) ;
  IntermediateNode() {_isSolved = false ;};
  
protected:
  
  bool _isSolved ;
  std::vector<TreeNode *> _children ; 
};