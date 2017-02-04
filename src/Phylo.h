#include <RcppArmadillo.h>
#include "GraphProperties.h"
#include "VertexProperties.h"
#include <boost/graph/adjacency_list.hpp>

using namespace Rcpp;
using namespace arma;
using namespace boost;

typedef adjacency_list<vecS, vecS, bidirectionalS, VertexProperties, no_property, GraphProperties, listS> DirectedGraph;
typedef graph_traits<DirectedGraph>::vertex_iterator vertex_iter;

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
  std::vector<cube> alignmentBin ;
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
  void writeTips() ;
  void updateInternalVertices() ;
  void updateVertex(const vertex_iter) ;
  graph_traits<DirectedGraph>::vertex_descriptor parentVertex(const vertex_iter v) ;
  // Output node children:
  uvec Children(const umat &, const uint) ;
  // The likelihood at one locus for one rate category.
  double logLikOneLocusOneRate(const uint, const int, const bool) ;
  vec getLogLikVec(const uint, mat &, const int) ;
  void internalFun(const uint, mat &, const int) ;
  void buildTreeGraph() ;
  void convertAlignmentList(const List &) ;
  
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