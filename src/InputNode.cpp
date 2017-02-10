#include <gsl/gsl_sf_pow_int.h>
#include "InputNode.h"

InputNode::InputNode(longVec alignmentBinOneTip, uint numRateCats)
{
  _solution.resize(numRateCats) ;
  
  for (auto & iter : _solution)
  {
    iter.resize(alignmentBinOneTip.size()) ;
    iter = alignmentBinOneTip ;
  }
}
   