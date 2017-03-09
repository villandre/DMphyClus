#include "InputNode.h"

// void InputNode::DeriveKey(solutionDictionaryType & solutionDictionary, const uint & rateCategory, const uint & matListIndex, const uint & locusNum)
// {
//   std::vector<uint> intEquiv = conv_to<std::vector<uint>>::from(_input->at(locusNum)) ;
//   std::vector<bool> boolEquiv(intEquiv.size()) ;
//   std::transform(intEquiv.begin(), intEquiv.end(), boolEquiv.begin(), [] (const uint integerValue) {return (bool) integerValue ;}) ;
//   _dictionaryKeyVec.at(locusNum) = std::hash<std::vector<bool>>{}(boolEquiv) ; // I should think collisions will not occur when hashing small vectors of booleans. After all, alignments, even with ambiguities, are strictly represented by vectors of 0's and 1's.
//   
//   if (solutionDictionary->at(rateCategory).count(_dictionaryKeyVec.at(locusNum)) == 0)
//   {
//     solutionDictionary->at(rateCategory)[_dictionaryKeyVec.at(locusNum)] = conv_to<vec>::from(*_input) ;
//   }
//   _keyDefined = true ;
// }

void InputNode::InitMapAndIterVec(solutionDictionaryType & solutionDictionary)
{
  std::vector<S> hashStructVec(_inputVec->size()) ; 
  std::transform(_inputVec->begin(), _inputVec->end(), hashStructVec.begin(), [] (const uvec & inputVec) 
  {
    std::vector<int> convertedVec = conv_to<std::vector<int>>::from(inputVec) ;
    std::vector<bool> convertedBool(convertedVec.size()) ;
    std::transform(convertedVec.begin(), convertedVec.end(), convertedBool.begin(), [] (int & i) {return (bool) i ;}) ;
    std::size_t hashedInput = std::hash<std::vector<bool>>{}(convertedBool) ;
    S myStruct = S(hashedInput, hashedInput, 0, 0) ;
    return myStruct ;
  }) ;
  std::pair<mapIterator, bool> insertResult ;
  for (unsigned int i = 0 ; i < (*_inputVec).size() ; i++)
  {
    for (uint j = 0 ; j < solutionDictionary->size(); j++)
    {
      std::pair<mapIterator, bool> insertResult = solutionDictionary->at(j).insert(std::pair<S,vec>(hashStructVec.at(i),conv_to<vec>::from(_inputVec->at(i)))) ;
    }
    _dictionaryIterVec.at(i) = insertResult.first ; 
  } // _dictionaryIterVec will point to elements in the last dictionary, but it doesn't matter, since the pointed solution is the input vector itself, which does not depend on the rate category
  
}