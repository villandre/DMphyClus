#include <gsl/gsl_sf_pow_int.h>
#include "InputNode.h"

void InputNode::DeriveKey(solutionDictionaryType & solutionDictionary = NULL) 
{
  std::vector<bool> boolEquiv = conv_to<std::vector<bool>>::from(_input) ; 
  _dictionaryKey = std::hash<std::vector<bool>>{}(boolEquiv) ; // I should think collisions will not occur when hashing small vectors of booleans. After all, alignments, even with ambiguities, are strictly represented by vectors of 0's and 1's.
}