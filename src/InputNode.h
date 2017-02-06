#include <boost/functional/hash.hpp>
#include <assert.h>
#include "TreeNode.h"

class InputNode:TreeNode
{
public:
  virtual bool IsSolved() {return true ;};
  virtual bool CanSolve() {return true ;};
  virtual void AddChild(TreeNode* child) {assert(false) ;};
  virtual void RemoveChild(TreeNode* child) {assert(false) ;};
  virtual void SetSolution(Col<long double> & inputVec) {_solution = inputVec ;};
  virtual void SetPattern() ;
  
protected:
  virtual void InvalidateSolution() {assert(false) ;};
};