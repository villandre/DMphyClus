// [[Rcpp::depends(RcppArmadillo)]]

#include <gperftools/profiler.h>
#include <RcppArmadillo.h>
#include <vector>
#include <array>
#include <stdexcept>
#include <omp.h>
#include <unordered_map>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <sys/resource.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

//#include <boost/algorithm/string.hpp>
// [[Rcpp::plugins(openmp)]]

using namespace arma;
using namespace Rcpp;

class phylo { // NOTE: ONCE CODE IS DEBUGGED, YOU CAN SPEED IT UP BY USING [] or .at() FOR INDEXING.
protected:
    umat edge ;
    uint Nnode ;
    uint numTips ;
    uint numLoci ;
    uint numUniqueLoci ;
    //umat alignment ;
    std::vector<mat> logTransMatVec ;
    vec logLimProbs ;
    uvec updateOrder ;
    double logLik ;
    std::vector<Rcpp::NumericMatrix> pruningMatVec ; // Each element of pruningMatVec comprises all the vectors at the root node used in computing the likelihood with the pruning algorithm. It follows that pruningMat[[x]] has one column per locus and one row per state (4 for nucleotide alignments, 20 for AA). This vector has as many elements as rate categories.
    //std::vector<mat> alignmentBinList ;
    List alignmentBin ;
    uint numStates;
    uint numRateCats;
    std::vector<std::vector<int> > children;
    int numOpenMPthreads ;
    //uvec childNodeInClusIndic ; // This is a vector with 0s and 1s. 1 means that the child node is at the end of an edge inside the cluster. 0 means that the child node is at the end of an edge outside a cluster. This is probably not functional and should be removed, until it's properly implemented.
    bool internalPhyloFlag ; // This flag indicates that this is an internal phylogeny, i.e. a phylogeny that supports subphylogenies for clusters. Concretely, this implies that the configuration corresponding to cluster centroids, that are at some of the tips of this phylogeny, are vectors with as many elements as rate categories.
    uvec sitePatterns ;
    // Member functions:
    // Compute the order in which nodes must be update using Felsenstein's tree-pruning algorithm.
    void compUpdateVec() ;
    // Output node children:
    uvec Children(const umat &, const uint) ;
    // Member function: The likelihood at one locus for one rate category.

    double logLikOneLocusOneRate(const uint, const int, const bool) ;
    vec getLogLikVec(const uint, mat &, const int) ;
    void internalFun(const uint, mat &, const int) ;

public:
    // constructor, sets up data structures
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
    //phylogeny() : numOpenMPthreads(4) {} ;
} ;

class phylogenyAlpha: public phylo {
public:
    // constructor, sets up data structures
    phylogenyAlpha(const NumericMatrix &, const CharacterMatrix &, const NumericVector &, const List &, const int, const CharacterVector &, const bool, const bool, const uvec &) ;
    phylogenyAlpha(const NumericMatrix &, const List &, const NumericVector &, const List &, const int, const bool, const bool, const uvec &) ;
    phylogenyAlpha() ; // Empty constructor
private:
    Rcpp::CharacterMatrix alignmentRcpp ;
    void convertSTL() ;
    void convertNucleoToNum() ;
    std::vector<std::vector<std::string>> alignmentAlpha ;
    std::vector<std::string> equivalency ; // This vector gives equivalencies between the nucleotides and numbers. This should match the rows and columns of the rate matrix, e.g. A = 1, T = 2, C = 3, G = 4;
    void defineMap() ;
    std::unordered_map<std::string, vec> mymap;
};

// The constructor

phylo::phylo(const NumericMatrix & edgeMat, const NumericVector & logLimProbsVec, const List & logTransMatList, const int numOpenMP, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) {

  mat edgeDouble = Rcpp::as<mat>(edgeMat); // The first argument (matrix) of the R function is called edgeMatrix.
  edge = conv_to<umat>::from(edgeDouble); // edge is a matrix of unsigned integers.
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
  
  children.resize(edge.n_rows + 2) ;
  
  for(uint i = 0; i < edge.n_rows; i++)
  {
    
    children[edge(i, 0)].push_back(edge(i, 1));
  }
  
  if (returnMatIndic) {
    
    pruningMatVec.resize(numRateCats) ;
    for (auto& i : pruningMatVec) {
      
      i = NumericMatrix(numStates, numLoci) ;
    }
  } else{}
}

phylogenyAlpha::phylogenyAlpha(const NumericMatrix & edgeMat, const CharacterMatrix & alignmentAlphaMat, const NumericVector & logLimProbsVec, const List & logTransMatList, const int numOpenMP, const CharacterVector & equivVec, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) : phylo(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsVec) {

    alignmentRcpp = alignmentAlphaMat ;
    equivalency = Rcpp::as<std::vector<std::string>>(equivVec) ;
    convertSTL() ;
    defineMap() ;
    convertNucleoToNum() ;
    // if (!(logLimProbsVec[0] == 1000.0)) { // For getConvertedAlignmentToWrao, edgeMat is a placeholder for a phylo with two tips, but alignmentAlphaMat.nrow() is not 2. This is why this last part won't run. I use a placeholder for logLimProbsVec too, with first element equal to 1000, hence this check.
    //   numTips = alignmentAlphaMat.nrow() ;
    //   Nnode = edgeMat.nrow() - numTips + 1 ;
    //   compUpdateVec() ;
    // } else{}
}

// This is the faster constructor, that does not convert the character alignment into matrices of identity vectors, because it's provided by the user (alignmentList).
phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::List & alignmentList, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numOpenMP, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) : phylo(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsVec) {
    
  alignmentBin = alignmentList ;
  if (!internalFlag) {
    
    numTips = as<NumericMatrix>(alignmentBin[0]).ncol() ;
    numUniqueLoci = alignmentBin.size() ;
  } else { // Tip vectors depend on the rate category.
    
    numTips = alignmentBin.size() ;
    numUniqueLoci = as<NumericMatrix>(as<List>(alignmentBin[0])[0]).ncol() ;
  }
  Nnode = edgeMat.nrow() - numTips + 1 ;
  compUpdateVec() ;
}

phylo::phylo() {}
phylogenyAlpha::phylogenyAlpha():phylo() {}

void phylogenyAlpha::convertSTL() {
    
    alignmentAlpha.resize(alignmentRcpp.ncol()) ;

    for (int i = 0; i < alignmentRcpp.ncol(); i++) {

        Rcpp::CharacterVector locusCol = alignmentRcpp(Rcpp::_, i) ;

        alignmentAlpha[i] = Rcpp::as<std::vector<std::string>>(locusCol) ;

    }
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
        if (myTest) {

            updateOrderVec(indexForUpdateVec) = numTips + i + 1 ;
            indexForUpdateVec += 1 ;
            readyForUpdate.push_back(i) ;
        } else{} ;
    }
    while (readyForUpdate[0] != 0) {

        uint parentOfReadyNode = parentVec(readyForUpdate[0]+numTips) ; //parentVec includes tips as well.
        uint parentOfReadyNodeIndex = parentOfReadyNode - numTips - 1;
        countForReady(parentOfReadyNodeIndex) += 1 ;
        bool myTest = (countForReady(parentOfReadyNodeIndex) == numChildren(parentOfReadyNodeIndex)) ;
        if (myTest) {

            readyForUpdate.push_back(parentOfReadyNodeIndex) ;
            updateOrderVec(indexForUpdateVec) = parentOfReadyNode ;
            indexForUpdateVec += 1 ;
        } else{} ;
        readyForUpdate.erase(readyForUpdate.begin()) ;
    }
    updateOrder = updateOrderVec ;
}

uvec phylo::Children(const umat & edgeMat, const uint parentNum) {

    uvec posVec = find(edgeMat.col(0) == parentNum) ;
    uvec secondColumn = edgeMat.col(1) ;
    return secondColumn.elem(posVec) ;
}

void phylogenyAlpha::defineMap() {

    for (uint i = 0; i < numStates; i++) {

        vec unitVec(numStates, fill::zeros) ;
        unitVec(i) = 1 ;
        std::string label = equivalency[i] ;
        mymap[label] = unitVec ;
    }
    if (numStates == 4) { // If we have four alphanumeric states, the assumption is that we're looking at a nucleotide alignment with states "a", "t", "c", "g".
    // Now, we handle ambiguities.
        mymap["r"] = (mymap["a"] + mymap["g"]) ;
        mymap["y"] = (mymap["c"] + mymap["t"]) ;
        mymap["s"] = (mymap["c"] + mymap["g"]) ;
        mymap["w"] = (mymap["a"] + mymap["t"]) ;
        mymap["k"] = (mymap["g"] + mymap["t"]) ;
        mymap["m"] = (mymap["a"] + mymap["c"]) ;
        mymap["b"] = (mymap["c"] + mymap["g"] + mymap["t"]) ;
        mymap["d"] = (mymap["a"] + mymap["g"] + mymap["t"]) ;
        mymap["h"] = (mymap["a"] + mymap["c"] + mymap["t"]) ;
        mymap["v"] = (mymap["a"] + mymap["c"] + mymap["g"]) ;
        mymap["-"] = (mymap["a"] + mymap["t"] + mymap["c"] + mymap["g"]) ;
        mymap["n"] = (mymap["a"] + mymap["t"] + mymap["c"] + mymap["g"]) ;
        mymap["."] = (mymap["a"] + mymap["t"] + mymap["c"] + mymap["g"]) ;
    } else{}
}

void phylogenyAlpha::convertNucleoToNum()  {
    
    List alignmentBinList(numLoci) ;
    std::vector<std::string> locusCol ;
    mat containerMat(numStates, alignmentRcpp.nrow(), fill::zeros) ;
    for (uint locusIndex = 0 ; locusIndex < numLoci ; locusIndex++) {
      
      locusCol = alignmentAlpha[locusIndex] ;
      for (uint i = 0; i < containerMat.n_cols; i++) {
        
        containerMat.col(i) = mymap[locusCol[i]] ;
      }
      alignmentBinList(locusIndex) = as<NumericMatrix>(wrap(containerMat)) ;
    }
    
    alignmentBin = alignmentBinList ; 
}

vec phylo::getLogLikVec(const uint childNodeIndex, mat & nodeTipMat, const int rateIndex) {

    double finiteMinLogLikAtTip ;
    rowvec logLikVec((nodeTipMat).n_rows) ;

    vec relocNodeTipVec(nodeTipMat.col(childNodeIndex - 1)); // Indices start at 0, not at 1.
    // We remove elements that are -Inf. They will equate 0 later on and so, are of no use.
    bool checkFinite = is_finite(relocNodeTipVec) ;
    vec finiteRelocNodeTipVec(relocNodeTipVec) ;
    mat transpFiniteTransMat(arma::trans(logTransMatVec[rateIndex])) ;
    if (!checkFinite) { // My guess would be that the subsetting is only required when childNodeIndex is inferior to numTips: in this case, the nodes are leaves.

        finiteRelocNodeTipVec = finiteRelocNodeTipVec.elem(find_finite(relocNodeTipVec)) ;
        transpFiniteTransMat = transpFiniteTransMat.rows(find_finite(relocNodeTipVec)) ;
    } else {}
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
  
  if (!internalPhyloFlag) {
    
    nodeTipMat.cols(0, numTips-1) = log(as<mat>(alignmentBin[locusNum])) ;
  } else {
    
    //nodeTipMat.cols(0, numTips-1) = multiAlignmentBinList[rateIndex][locusNum] ; // The alignment now depends on the rate category.
    
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
  
  if (returnMatIndic) {
    
    for (uint i = 0 ; i < nodeTipMat.n_rows; i++) {
      double valueInMatrix = nodeTipMat.at(i,numTips) ;
      pruningMatVec[rateIndex](i, locusNum) = valueInMatrix ;
    }
  } else{}
  
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
  if (!returnMatIndic) {
    vec rowMin = min(logLikMat,1) ;
    logLikMat.each_col() -= rowMin ;
    vec logLiksToSum = rowMin + log(sum(exp(logLikMat), 1)) - log(numRateCats) ;
    logLik = sum(logLiksToSum.elem(sitePatterns - 1)) ; // Indices begin at zero hence the -1...
    //logLik = minLogLikSumRate + sum(exp(logLikSumRate - minLogLikSumRate)) - log(numRateCats)  ;
  } else{}
}

template<class Mat>
void print_matrix(Mat matrix) {

    matrix.print(std::cout);
}


template<class Col>
void print_vector(Col colvec) {

    colvec.print(std::cout);
}

//provide explicit instantiations of the template function for
//every matrix type you use somewhere in your program.
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);
template void print_vector<arma::uvec>(arma::uvec colvec);
template void print_vector<arma::vec>(arma::vec colvec);

// [[Rcpp::export]]

SEXP logLikCppToWrap(NumericMatrix & edgeMat, NumericVector & logLimProbsVec, List & logTransMatList, int numOpenMP, SEXP & equivVector, List alignmentBin, const bool returnMatIndic, const bool internalFlag, const NumericVector sitePatterns) {

    omp_set_num_threads(numOpenMP) ;
    //uvec sitePatternsAmended(as<List>(alignmentBin[0]).size()) ;
    uint numLoci ;
    uvec sitePatternsAmended ;
    if (internalFlag) { 
      numLoci = as<NumericMatrix>(as<List>(alignmentBin[0])[0]).ncol() ;
      sitePatternsAmended.resize(numLoci) ;
      sitePatternsAmended = conv_to<uvec>::from(linspace(1, numLoci, numLoci))  ;
    } else{

      sitePatternsAmended = as<uvec>(sitePatterns) ;
    }
    SEXP container ;
   
    //ProfilerStart("/home/villandre/profileOut.out") ;
    phylogenyAlpha phyloObject(edgeMat, alignmentBin, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsAmended) ;
    phyloObject.logLikPhylo(returnMatIndic);
    //ProfilerStop() ;
    if (returnMatIndic) {

      container = phyloObject.getPruningMat();
    } else {

      container = phyloObject.getLogLik();

    }
    return container ;
}

// [[Rcpp::export]]

SEXP getConvertedAlignmentToWrap(int numOpenMP, SEXP & equivVector, Rcpp::CharacterMatrix & alignmentAlphaMat, NumericVector & sitePatterns) {
  
  double src[] = {3, 3, 1, 2};
  NumericMatrix placeholderMat(2, 2, src) ; //src is considered an iterator.
  NumericVector placeholderVec(as<CharacterVector>(equivVector).size(),1000.0) ;
  List placeholderList(1) ; // This is a placeholder: 1 is an arbitrary number.
  placeholderList[0] = Rcpp::NumericMatrix(2,2) ;
  bool returnMatIndic = false ;
  bool internalFlag = false ;
  
  //phylogenyAlpha phyloObject = phylogenyAlpha(placeholderMat, alignmentAlphaMat, placeholderVec, placeholderList, numOpenMP, equivVector, placeholderVec, returnMatIndic, internalFlag, sitePatterns) ;
  phylogenyAlpha phyloObject(placeholderMat, alignmentAlphaMat, placeholderVec, placeholderList, numOpenMP, equivVector, returnMatIndic, internalFlag, as<uvec>(sitePatterns)) ;
  return phyloObject.getAlignmentBin() ;
}

template<class... Conts>
auto zip_range(Conts&... conts)
  -> decltype(boost::make_iterator_range(
      boost::make_zip_iterator(boost::make_tuple(conts.begin()...)),
      boost::make_zip_iterator(boost::make_tuple(conts.end()...))))
  {
    return {boost::make_zip_iterator(boost::make_tuple(conts.begin()...)),
            boost::make_zip_iterator(boost::make_tuple(conts.end()...))};
  }

// [[Rcpp::export]]

List getSitePatterns(List alignmentBin) {
  
  std::size_t myHash ;
  std::vector<std::size_t> hashVector ;
  NumericVector sitePatterns(alignmentBin.size()) ;
  std::vector<NumericMatrix> uniqueDNAdataBin ;
  auto hashIterator = hashVector.begin() ;
  std::vector<NumericMatrix> convertedAlignmentBin = as<std::vector<NumericMatrix> >(alignmentBin) ;
  
  for(auto&& i : zip_range(convertedAlignmentBin, sitePatterns)) {
    
    myHash = boost::hash_range(as<NumericMatrix>(i.get<0>()).begin(), as<NumericMatrix>(i.get<0>()).end()) ;
    
    hashIterator = std::find(hashVector.begin(), hashVector.end(), myHash) ;
    i.get<1>() = std::distance(hashVector.begin(), hashIterator) + 1 ;

    if (hashIterator == hashVector.end()) {

      hashVector.push_back(myHash) ;
      uniqueDNAdataBin.push_back(as<NumericMatrix>(i.get<0>())) ;
    }
  }
  
  return List::create(Named("uniqueDNAdataBin") = wrap(uniqueDNAdataBin),
                      Named("sitePatterns") = wrap(sitePatterns));
}

// [[Rcpp::export]]aa

// double clusIndLogPrior(Rcpp::NumericVector clusInd, double alpha) { // Only works when there are no gaps in clusInd.
//
//     std::vector<int> clusCounts(Rcpp::max(clusInd), 0) ;
//     int numElements = clusInd.size() ;
//
//     clusCounts.at(clusInd[0]-1) = 1 ;
//
//     arma::vec logContributions(numElements, fill::zeros) ;
//
//     double probNumer ;
//     for (int i = 1; i < numElements; i++) { // It must begin at 1 instead of 0 because element 0 contributes 1 trivially.
//        probNumer = clusCounts.at(clusInd.at(i)-1) ;
//         if (probNumer == 0) {
//             probNumer = alpha ;
//         } else{}
//         clusCounts.at(clusInd.at(i)-1) = clusCounts.at(clusInd.at(i) - 1) + 1 ;
//         logContributions.at(i) = log(probNumer) - log(i+alpha) ;
//     }
//
//     return(sum(logContributions)) ;
// }

// // [[Rcpp::exportaaa]]
//
// double clusIndLogPrior(Rcpp::NumericVector clusInd, double alpha) { // Only works when there are no gaps in clusInd.
//
//     std::vector<int> clusCounts(Rcpp::max(clusInd), 0) ;
//     int numElements = clusInd.size() ;
//     int numClusters = 1 ; //The trivial number of clusters in the vector made up only of the first element of clusInd.
//
//     clusCounts.at(clusInd[0]-1) = 1 ;
//
//     arma::vec logContributions(numElements, fill::zeros) ;
//
//     double probNumer ;
//     int myCount ;
//     for (int i = 1; i < numElements; i++) { // It must begin at 1 instead of 0 because element 0 contributes 1 trivially.
//        myCount = clusCounts.at(clusInd.at(i)-1) ;
//        probNumer = 1 ;
//         if (myCount == 0) {
//             probNumer = alpha ;
//             numClusters = numClusters + 1 ;
//         } else{}
//         clusCounts.at(clusInd.at(i)-1) = clusCounts.at(clusInd.at(i) - 1) + 1 ;
//         logContributions.at(i) = log(probNumer) - log(numClusters+alpha) ;
//     }
//
//     return(sum(logContributions)) ;
// }
