#include "InputNode.h"

void InputNode::DeriveKey(solutionDictionaryType & solutionDictionary)
{
  std::vector<uint> intEquiv = conv_to<std::vector<uint>>::from(_input) ;
  std::vector<bool> boolEquiv(intEquiv.size()) ;
  std::transform(intEquiv.begin(), intEquiv.end(), boolEquiv.begin(), [] (const uint integerValue) {return (bool) integerValue ;}) ;
  _dictionaryKey = std::hash<std::vector<bool>>{}(boolEquiv) ; // I should think collisions will not occur when hashing small vectors of booleans. After all, alignments, even with ambiguities, are strictly represented by vectors of 0's and 1's.
  _keyDefined = true ;
}