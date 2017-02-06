#include <gsl/gsl_sf_pow_int.h>
#include "InputNode.h"

void InputNode::SetPattern() 
{
  uint myHash = 0;
  for (int i = 0; i < _solution.size(); i++) {
    myHash = myHash + gsl_sf_pow_int(10, i)*_solution.at(i) ;
  }
  _pattern = myHash ;
}
   