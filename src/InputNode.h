#include "TreeNode.h"

class InputNode:public TreeNode
{
public:
  
  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChild(TreeNode *) {assert(false) ;}
  std::vector<TreeNode *> GetChildren() {std::vector<TreeNode *> myVec; myVec.push_back(NULL) ; return myVec;} // An input node returns a null pointer when it is asked to provide the address of a child.
  void RemoveChildren() {}
  
  bool IsSolved() {return true ;}
  bool CanSolve() {return true ;}
  void SetSolved(bool status) {}
  void ComputeSolutions(solutionDictionaryType &, const std::vector<mat> &, const uint &, boost::asio::io_service *, boost::mutex &, boost::barrier &) {assert(false) ;}
  bool ComputeSolution(solutionDictionaryType & dictionary, const std::vector<mat> &, const uint &, const uint &, boost::mutex &) {assert(false) ; return true ;} //Solution is known, this should not get called.
  void InvalidateSolution() {assert(false) ;}
  
  void SetInput(std::vector<uvec> * inputVec) { _inputVec = inputVec ;}
  std::vector<uvec> * GetInput() { return _inputVec ;}
  
  void InitMapAndIterVec(solutionDictionaryType &, const uint &) ;
  void RestoreIterVecAndExp() {} // Solutions for input nodes are trivial and never change. It follows that a restore should not do anything.
  mapIterator GetDictionaryIterator(const uint & elementNum, const uint & numRateCats) {return _dictionaryIterVec.at(elementNum) ;}
  S GetSfromVertex(const uint &, const uint &, const uint &) {assert(false) ;};
  
  InputNode(uint & numLoci, uint & numRates, solutionDictionaryType & solutionDictionary)
  {
    _parent = NULL ;
    _dictionaryIterVec.resize(numLoci) ;
    for (uint i = 0 ; i < numLoci ; i++)
    {
      _dictionaryIterVec.at(i) = solutionDictionary->begin() ;
    }
    _updateFlag = false ;
  }
  
protected:
  std::vector<uvec> * _inputVec ;
};
