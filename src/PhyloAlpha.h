#include "Phylo.h"
using namespace Rcpp;
using namespace arma;

class phylogenyAlpha: public phylo {

  // Attributes
private:

  std::unordered_map<std::string, vec> mymap;
  CharacterMatrix alignmentRcpp ;
  std::vector<std::vector<std::string>> alignmentAlpha ;
  std::vector<std::string> equivalency ; // This vector gives equivalencies between the nucleotides and numbers. This should match the rows and columns of the rate matrix, e.g. A = 1, T = 2, C = 3, G = 4;
  
// Methods

private:
  void convertSTL() ;
  void convertNucleoToNum() ;
  void defineMap() ;
  void initializeGraph() ;
  
public:
  phylogenyAlpha(const NumericMatrix &, const CharacterMatrix &, const NumericVector &, const List &, const int, const CharacterVector &, const bool, const bool, const uvec &) ;
  phylogenyAlpha(const NumericMatrix &, const List &, const NumericVector &, const List &, const int, const bool, const bool, const uvec &) ;
  phylogenyAlpha() ; // Empty constructor
};
