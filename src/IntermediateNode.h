#include <boost/functional/hash.hpp>
#include <assert.h>
#include "TreeNode.h"

class mapKey {
public:
  TreeNode * _currentNode ;
  bool _nodeWithinClus ;
  uint _transProbIndexAtNode ;
  uint _rateCategory ;
  
  bool operator==(const mapKey& myMap) const
  {
    return true ;
  } ;
  
  mapKey(TreeNode * currentNode, bool nodeWithinClus, uint transProbIndex): _currentNode(currentNode), _nodeWithinClus(nodeWithinClus), _transProbIndexAtNode(transProbIndex) {} ;
};

// custom hash
struct MyHash
{
  std::size_t operator()(mapKey const& key) const 
  {
    std::size_t aHash = 2 ;
    // aHash = aHash^key._child1Key ;
    // aHash = aHash^key._child2Key ;
    // aHash+= (std::size_t) key._nodeWithinClus ;
    return aHash ; 
  }
};


class IntermediateNode : public TreeNode
{
public:
  
  bool IsSolved() {return _isSolved ;};
  bool CanSolve() ;
  std::size_t GetPattern() ; // size_t * are pointers to the hash keys attributed to the children nodes by boost::hash. Beware of collisions!
  void AddChild(TreeNode*) ;
  void RemoveChild(TreeNode*) ;
  void SetSolution(Col<long double> &) ;
  void InvalidateSolution() ;
  void SetPattern(std::unordered_map<mapKey, Col<long double>, MyHash> &) ;
  void ToggleSolved() {_isSolved = !_isSolved ;};
  
  IntermediateNode(): _isSolved(false)  {_parent = NULL ;};
  
protected:
  
  bool _isSolved ;
  std::vector<TreeNode *> _children ;
};