#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

struct VertexProperties 
{ 
  unsigned id;
  Cube<double> CubeLiOneSlicePerRate;
  bool valueAssigned = FALSE ;
  VertexProperties() : id(0) {}
  VertexProperties(unsigned i) : id(i) {}
  uint numChildrenDefined = 0 ; 
};
