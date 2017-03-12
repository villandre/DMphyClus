#include "InputNode.h"

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
  for (unsigned int i = 0 ; i < _inputVec->size() ; i++)
  {
    for (uint j = 0 ; j < solutionDictionary->size(); j++)
    {
      insertResult = solutionDictionary->at(j).insert(std::pair<S,std::pair<vec,float>>(hashStructVec.at(i),std::pair<vec,float>(conv_to<vec>::from(_inputVec->at(i)), 0))) ;
      _dictionaryIterVec.push_back(insertResult.first) ;
    }
  }
}