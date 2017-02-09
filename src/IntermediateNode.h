#include <boost/functional/hash.hpp>
#include <assert.h>
#include "TreeNode.h"

class mapKey {
public:
  TreeNode * _currentNode ;
  bool _nodeWithinClus ;
  uint _transProbIndexAtNode ;
  
  bool operator==(const mapKey& myMap) const
  {
    return (this->_nodeWithinClus == myMap._nodeWithinClus) && 
      (this->_transProbIndexAtNode == myMap._transProbIndexAtNode) &&
      (((this->_currentNode->GetChildren()[0]->GetPattern() == myMap._currentNode->GetChildren()[0]->GetPattern()) && 
      (this->_currentNode->GetChildren()[1]->GetPattern() == myMap._currentNode->GetChildren()[1]->GetPattern())) || 
      ((this->_currentNode->GetChildren()[0]->GetPattern() == myMap._currentNode->GetChildren()[1]->GetPattern()) && 
      (this->_currentNode->GetChildren()[1]->GetPattern() == myMap._currentNode->GetChildren()[0]->GetPattern()))) ;
  } ;
  
  mapKey(TreeNode * currentNode, bool nodeWithinClus, uint transProbIndex): _currentNode(currentNode), _nodeWithinClus(nodeWithinClus), _transProbIndexAtNode(transProbIndex) {} ;
};

// custom hash
struct MyHash
{
  std::size_t operator()(mapKey const& key) const 
  {
    std::size_t aHash = 2 ;
    aHash = aHash^key._child1Key ;
    aHash = aHash^key._child2Key ;
    aHash+= (std::size_t) key._nodeWithinClus ;
    return aHash; 
  }
};


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
  virtual void SetPattern(std::unordered_map<mapKey, Col<long double>, MyHash> &) ;
  IntermediateNode() {_isSolved = false ; _parent = NULL ;};
  virtual void ToggleSolved() {_isSolved = !_isSolved ;};
  
protected:
  
  bool _isSolved ;
  std::vector<TreeNode *> _children ; 
};