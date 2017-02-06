//#include <RcppArmadillo.h>
#include "Phylo.h"
using namespace Rcpp;
using namespace arma;
using namespace boost;

phylo::phylo(const NumericMatrix & edgeMat, const NumericVector & logLimProbsVec, const List & logTransMatList, const int numOpenMP, const bool returnMatIndic, const uvec & sitePatternsVec) {
  
  //mat edgeDouble = Rcpp::as<mat>(edgeMat); // The first argument (matrix) of the R function is called edgeMatrix.
  edge = as<umat>(edgeMat) ;
  
  //edge = conv_to<umat>::from(edgeDouble); // edge is a matrix of unsigned integers.
  logLimProbs = Rcpp::as<vec>(logLimProbsVec); // The second argument (numeric) of the R function is called logLimProbs.
  numRateCats = logTransMatList.size() ;
  numStates = logLimProbsVec.size() ;
  numOpenMPthreads = numOpenMP ;
  internalPhyloFlag = internalFlag ;
  sitePatterns = sitePatternsVec ;
  numLoci = sitePatterns.n_rows ;
  numUniqueLoci = sitePatterns.max() ;
  
  logTransMatVec.resize(numRateCats) ;
  
  std::transform(logTransMatList.begin(), logTransMatList.end(), logTransMatVec.begin(), [] (const NumericMatrix &initialMatrix) {return as<mat>(initialMatrix);}) ;
  
  
  if (returnMatIndic) 
  {
    pruningMatVec.resize(numRateCats) ;
    for (auto& i : pruningMatVec) {
      
      i = NumericMatrix(numStates, numLoci) ;
    }
  } 
}

// void phylo::convertAlignmentList(const List & alignmentList) {
//   std::vector<NumericVector> foo=as<std::vector<NumericVector>>(alignmentList) ;
//   cube cubeArray(as<NumericVector>(alignmentList[0]).begin(), 4, 4, 4, false);
//   std::transform(foo.begin(), foo.end(), alignmentBin.begin(), [](NumericVector & a) { // Why can't I have a const here?
//     IntegerVector arrayDims = a.attr("dim");
//     cube cubeArray(a.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
//     return conv_to<Cube<long double> >::from(cubeArray) ;
//   }) ;
// }

phylo::phylo() {}

void phylo::buildTreeGraph() {
  
  typedef graph_traits<DirectedGraph>::vertices_size_type Vsize;
  typedef std::vector<std::pair<uint,uint>> pairVector ;
  edge = edge - 1 ;// Vertices must be encoded from 0 to n.
  pairVector edgeRecast(edge.n_rows) ;
  std::transform(edge.begin(), edge.begin() + edge.n_rows, edge.begin() + edge.n_rows + 1, edgeRecast.begin(),
                 [](const int& a, const int& b) { return std::make_pair(a, b); });
  DirectedGraph myGraph(edgeRecast.begin(), edgeRecast.end(), (Vsize) edge.max()) ;
  phyloGraph = myGraph ;
}

graph_traits<DirectedGraph>::vertex_descriptor phylo::parentVertex(const vertex_iter v) {
  return source(*in_edges(*v, phyloGraph).first, phyloGraph) ;
} 

void phylo::initializeVertices() {
  
  // get the property map for vertex indices
  typedef property_map<DirectedGraph, vertex_index_t>::type IndexMap;
  IndexMap index = get(vertex_index, phyloGraph);
  
  std::pair<vertex_iter, vertex_iter> vp = vertices(phyloGraph);
  
  // Takes advantage of the fact that tips are identified by numbers 0, ..., numTips - 1.
  
  for (vertex_iter v = vp.first; v < vp.first + alignmentBin.size();  v++) { // Since lambda functions are basically closures, I don't see how I can use for_each in this setting.
    phyloGraph[*v].partLikVecOneElementPerRate = alignmentBin[index[*v]] ;
    phyloGraph[*v].valueAssigned = TRUE ;
    phyloGraph[parentVertex(v)].numChildrenDefined++ ;
  }
}


void phylo::updateInternalVertices() {
  
  std::vector<vertex_iter> updateVertexIterVec ;
  updateVertexIterVec.reserve((int) num_vertices(phyloGraph) - numTips) ; // This is to prevent re-allocation when the vector is grown, which would mess up the iterator.
  
  // We find a place to start the update by doing a first pass through the internal nodes.
  
  for (vertex_iter v = vertices(phyloGraph).first + numTips ; v != vertices(phyloGraph).second; v++) {
    if (phyloGraph[*v].numChildrenDefined >= out_degree(*v, phyloGraph)) { // Since internal vertices values will be saved and restored later on, it might be that numChildrenDefined ends up larged than the number of children.
      updateVertexIterVec.push_back(v) ;
    }
  }
  
  // Now, we actually perform the updates.
   
  std::vector<vertex_iter>::iterator currentPos = updateVertexIterVec.begin() ; 
  while (currentPos != updateVertexIterVec.end()) { // Hopefully, .end() does not refer to the capacity, but to the actual number of elements in a vector...
    if (!phyloGraph[**currentPos].valueAssigned) {
      updateVertex(*currentPos) ;
      phyloGraph[**currentPos].valueAssigned = TRUE ;
    }
    if (in_degree(**currentPos, phyloGraph) == 0) { // We have reached the root, we must break the loop.
      break ;
    } // The while condition is useless... This logical scheme should be reformulated.
    phyloGraph[parentVertex(*currentPos)].numChildrenDefined++ ;
    if (phyloGraph[parentVertex(*currentPos)].numChildrenDefined >= out_degree(parentVertex(*currentPos), phyloGraph)) {
      updateVertexIterVec.push_back(*currentPos) ;
    }
    currentPos++ ;// currentPos indicates the position in the updateVertexIterVec. Once no new element can be added, it will move to .end(), thus ending the loop.
  }
}

// IMPORTANT FUNCTION: THIS IS AN INTEGRAL PART OF THE TREE_PRUNING ALGORITHM
void phylo::updateVertex(const vertex_iter v) {
  
  phyloGraph[*v].partLikVecOneElementPerRate = 
}

void phylo::compUpdateVec() {
  
  uint numNodesTotal = edge.n_rows + 1 ;
  uvec updateOrderVec(Nnode, fill::zeros) ;
  uvec numChildren(Nnode) ;
  
  auto childrenIter = std::next(children.begin(), numTips + 1) ;
  std::transform(childrenIter, std::next(childrenIter, Nnode), numChildren.begin(), [] (std::vector<int> &arg) { return arg.size(); }) ;
  
  uvec parentVec(numNodesTotal, fill::zeros) ;
  for (uint i = 0; (i < edge.n_rows); i ++) {
    
    parentVec(edge(i,1) - 1) = edge(i,0) ;
  } // The root will have a 0, indicating it has no parent.
  uvec countForReady(Nnode, fill::zeros) ;
  std::vector<uint> readyForUpdate ;
  for (uint i = 0; (i < numTips); i++) {
    
    countForReady(parentVec(i) - numTips - 1) += 1 ;   //countForReady is only for internalNodes, hence the subtraction of (numTips + 1)
  }
  uint indexForUpdateVec = 0 ;
  for (uint i = 0; (i < Nnode) ; i++) {
    
    bool myTest = (countForReady(i) == numChildren(i)) ;
    if (myTest) 
    {
      updateOrderVec(indexForUpdateVec) = numTips + i + 1 ;
      indexForUpdateVec += 1 ;
      readyForUpdate.push_back(i) ;
    }
  }
  while (readyForUpdate[0] != 0) {
    
    uint parentOfReadyNode = parentVec(readyForUpdate[0]+numTips) ; //parentVec includes tips as well.
    uint parentOfReadyNodeIndex = parentOfReadyNode - numTips - 1;
    countForReady(parentOfReadyNodeIndex) += 1 ;
    bool myTest = (countForReady(parentOfReadyNodeIndex) == numChildren(parentOfReadyNodeIndex)) ;
    if (myTest) 
    {
      readyForUpdate.push_back(parentOfReadyNodeIndex) ;
      updateOrderVec(indexForUpdateVec) = parentOfReadyNode ;
      indexForUpdateVec += 1 ;
    }
    readyForUpdate.erase(readyForUpdate.begin()) ;
  }
  updateOrder = updateOrderVec ;
}

uvec phylo::Children(const umat & edgeMat, const uint parentNum) {
  
  uvec posVec = find(edgeMat.col(0) == parentNum) ;
  uvec secondColumn = edgeMat.col(1) ;
  return secondColumn.elem(posVec) ;
}

vec phylo::getLogLikVec(const uint childNodeIndex, mat & nodeTipMat, const int rateIndex) {
  
  double finiteMinLogLikAtTip ;
  rowvec logLikVec((nodeTipMat).n_rows) ;
  
  vec relocNodeTipVec(nodeTipMat.col(childNodeIndex - 1)); // Indices start at 0, not at 1.
  // We remove elements that are -Inf. They will equate 0 later on and so, are of no use.
  bool checkFinite = is_finite(relocNodeTipVec) ;
  vec finiteRelocNodeTipVec(relocNodeTipVec) ;
  mat transpFiniteTransMat(arma::trans(logTransMatVec[rateIndex])) ;
  if (!checkFinite) 
  { // My guess would be that the subsetting is only required when childNodeIndex is inferior to numTips: in this case, the nodes are leaves.
    finiteRelocNodeTipVec = finiteRelocNodeTipVec.elem(find_finite(relocNodeTipVec)) ;
    transpFiniteTransMat = transpFiniteTransMat.rows(find_finite(relocNodeTipVec)) ;
  } 
  finiteMinLogLikAtTip = min(finiteRelocNodeTipVec) ;
  mat repRelocVecMat(finiteRelocNodeTipVec.size(), relocNodeTipVec.size()) ;
  repRelocVecMat.each_col() = finiteRelocNodeTipVec ;
  logLikVec = finiteMinLogLikAtTip + log(sum(exp(transpFiniteTransMat + repRelocVecMat - finiteMinLogLikAtTip), 0)) ; // If logTransMatVec is always used in its transposed form, better do the transposition only once, in the R function.
  
  return arma::conv_to<vec>::from(logLikVec) ; // A column vector...
}

void phylo::internalFun(const uint parentNum, mat & nodeTipMat, const int rateIndex) {
  
  vec logLikVecByChild(nodeTipMat.n_rows, fill::zeros) ;\
  uint numChildren = children[parentNum].size() ;
  uvec theChildren(numChildren) ;
  for(uint i = 0; i < numChildren; i++) {
    
    theChildren(i) = children[parentNum][i];
  }
  for (uint i = 0 ; (i < numChildren); i++) {
    
    logLikVecByChild = logLikVecByChild + getLogLikVec(theChildren(i), nodeTipMat, rateIndex) ;
  }
  
  nodeTipMat.col(parentNum - 1) = logLikVecByChild ;
}

double phylo::logLikOneLocusOneRate(const uint locusNum, const int rateIndex, const bool returnMatIndic) {
  
  mat nodeTipMat = zeros<mat>(numStates, edge.max()) ; // This is a matrix that stores likelihood vectors. Each vector has as many elements as potential states.
  
  if (!internalPhyloFlag) 
  {
    nodeTipMat.cols(0, numTips-1) = log(as<mat>(alignmentBin[locusNum])) ;
  } 
  else
  {
    for (uint i = 0 ; i < numTips; i++) {
      
      NumericVector theColumn = as<NumericMatrix>(as<List>(alignmentBin[i])[rateIndex])(_, locusNum) ;
      nodeTipMat.col(i) = as<vec>(theColumn) ;
    }
  }
  for (uint i = 0; (i < updateOrder.n_rows) ; i++) {
    
    internalFun(updateOrder(i), nodeTipMat, rateIndex) ;
  }
  
  double minLogLikRoot = min(nodeTipMat.col(numTips)) ;
  
  double logLikForOneLocusOneRate = minLogLikRoot + log(sum(exp(nodeTipMat.col(numTips) + logLimProbs - minLogLikRoot))) ;
  
  if (returnMatIndic) 
  {
    for (uint i = 0 ; i < nodeTipMat.n_rows; i++) {
      double valueInMatrix = nodeTipMat.at(i,numTips) ;
      pruningMatVec[rateIndex](i, locusNum) = valueInMatrix ;
    }
  }
  
  return logLikForOneLocusOneRate ;
}

SEXP phylo::getAlignmentBin() {
  
  return(wrap(alignmentBin)) ;
}

// umat phylo::getAlignment() {
//     return alignment ;
// }

uint phylo::getNnode() {
  
  return Nnode ;
}

uint phylo::getNumStates() {
  
  return numStates ;
}

uint phylo::getNumRateCats() {
  
  return numRateCats ;
}

SEXP phylo::getLogLik() {
  
  return Rcpp::wrap(logLik) ;
}

uint phylo::getNumTips() {
  
  return numTips ;
}

SEXP phylo::getPruningMat() {
  
  std::vector<NumericMatrix> pruningMatVecReadj(pruningMatVec.size()) ;
  for (uint i = 0 ; i < pruningMatVec.size(); i++) {
    
    pruningMatVecReadj[i] = as<NumericMatrix>(wrap(conv_to<mat>::from(as<mat>(pruningMatVec[i]).cols(sitePatterns - 1)))) ; 
  }
  return wrap(pruningMatVecReadj) ;
}

void phylo::logLikPhylo(const bool returnMatIndic) {
  
  mat logLikMat(numUniqueLoci, numRateCats, fill::ones);
  
  for (uint rateNum = 0; (rateNum < numRateCats); rateNum++) {
    
    //#pragma omp parallel for
    for(uint locusNum = 0; locusNum < numUniqueLoci; locusNum++) {
      
      logLikMat.at(locusNum, rateNum) = logLikOneLocusOneRate(locusNum, rateNum, returnMatIndic) ;
    }
  }
  if (!returnMatIndic) 
  {
    vec rowMax = max(logLikMat,1) ;
    logLikMat.each_col() -= rowMax ;
    vec logLiksToSum = rowMax + log(sum(exp(logLikMat), 1)) - log(numRateCats) ;
    logLik = sum(logLiksToSum.elem(sitePatterns - 1)) ; // Indices begin at zero hence the -1...
  }
}
