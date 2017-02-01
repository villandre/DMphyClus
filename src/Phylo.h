#include <RcppArmadillo.h>
#include <boost/graph/adjacency_list.hpp>
using namespace Rcpp;
using namespace arma;

struct VertexProperties 
{ 
  unsigned id;
  Cube<double> intermediateCube;
  VertexProperties() : id(0) {}
  VertexProperties(unsigned i) : id(i) {}
};

struct GraphProperties 
{ 
  unsigned id ;
  GraphProperties() ; // Empty constructor
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::edge_iterator edge_iterator;

class phylo { // NOTE: ONCE CODE IS DEBUGGED, YOU CAN SPEED IT UP BY USING [] or .at() FOR INDEXING.

  // Attributes
  
protected:
  DirectedGraph phyloGraph ;
  umat edge ;
  uint Nnode ;
  uint numTips ;
  uint numLoci ;
  uint numUniqueLoci ;
  std::vector<mat> logTransMatVec ;
  vec logLimProbs ;
  uvec updateOrder ;
  double logLik ;
  std::vector<NumericMatrix> pruningMatVec ; // Each element of pruningMatVec comprises all the vectors at the root node used in computing the likelihood with the pruning algorithm. It follows that pruningMat[[x]] has one column per locus and one row per state (4 for nucleotide alignments, 20 for AA). This vector has as many elements as rate categories.
  List alignmentBin ;
  uint numStates;
  uint numRateCats;
  std::vector<std::vector<int> > children;
  int numOpenMPthreads ;
  bool internalPhyloFlag ; // This flag indicates that this is an internal phylogeny, i.e. a phylogeny that supports subphylogenies for clusters. Concretely, this implies that the configuration corresponding to cluster centroids, that are at some of the tips of this phylogeny, are vectors with as many elements as rate categories.
  uvec sitePatterns ;
  
  // Methods
  
protected:
  // Compute the order in which nodes must be update using Felsenstein's tree-pruning algorithm.
  void compUpdateVec() ;
  // Output node children:
  uvec Children(const umat &, const uint) ;
  // The likelihood at one locus for one rate category.
  double logLikOneLocusOneRate(const uint, const int, const bool) ;
  vec getLogLikVec(const uint, mat &, const int) ;
  void internalFun(const uint, mat &, const int) ;
  void buildTreeGraph() ;
  
public:
  phylo(const NumericMatrix &, const NumericVector &, const List &, const int, const bool, const bool, const uvec &) ;
  phylo() ;
  void logLikPhylo(const bool) ;
  SEXP getAlignmentBin() ;
  umat getAlignment() ;
  uint getNnode() ;
  uint getNumStates() ;
  uint getNumRateCats() ;
  uint getNumTips() ;
  SEXP getLogLik() ;
  SEXP getPruningMat() ;
} ;