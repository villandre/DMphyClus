#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

struct VertexProperties 
{ 
  uint _id;
  bool _isSolved = false;
  
  
  std::vector<Mat<long double>> partLikVecOneElementPerRate ;
  bool valueAssigned = FALSE ;
  VertexProperties() : _id(0) {}
  VertexProperties(unsigned i) : _id(i) {}
  uint numChildrenDefined = 0 ; 
};
