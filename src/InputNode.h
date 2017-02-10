#include <boost/functional/hash.hpp>
#include <assert.h>
#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  bool IsSolved() {return true ;};
  bool CanSolve() {return true ;};
  void AddChild(TreeNode* child) {assert(false) ;};
  void RemoveChild(TreeNode* child) {assert(false) ;};
  void SetSolution(Col<long double> & inputVec) {assert(false) ;};
  void GetSolution() ;
  void ToggleSolved() {};
  void SetPattern() ;
  
  InputNode(longVec, uint) ; // An input node has inputs structured as a vector, with each element corresponding to a different locus.
  
protected:
  void InvalidateSolution() {assert(false) ;};
  std::vector<std::vector<uint>> _input ;
};