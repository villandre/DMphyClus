#include "InputNode.h"

void InputNode::DeriveKey(solutionDictionaryType & solutionDictionary, const uint & rateCategory, const uint & matListIndex, const uint & locusNum)
{
  std::vector<uint> intEquiv = conv_to<std::vector<uint>>::from(_input->at(locusNum)) ;
  std::vector<bool> boolEquiv(intEquiv.size()) ;
  std::transform(intEquiv.begin(), intEquiv.end(), boolEquiv.begin(), [] (const uint integerValue) {return (bool) integerValue ;}) ;
  _dictionaryKeyVec.at(locusNum) = std::hash<std::vector<bool>>{}(boolEquiv) ; // I should think collisions will not occur when hashing small vectors of booleans. After all, alignments, even with ambiguities, are strictly represented by vectors of 0's and 1's.
  
  if (solutionDictionary->at(rateCategory).count(_dictionaryKeyVec.at(locusNum)) == 0)
  {
    solutionDictionary->at(rateCategory)[_dictionaryKeyVec.at(locusNum)] = conv_to<vec>::from(*_input) ;
  }
  _keyDefined = true ;
}
