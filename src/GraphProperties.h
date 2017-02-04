#include <RcppArmadillo.h>
using namespace arma;

struct GraphProperties 
{ 
  cube transProbMatCube ;
  GraphProperties() ; // Empty constructor
  GraphProperties(mat) ; 
};
