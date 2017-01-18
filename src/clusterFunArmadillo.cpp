// [[Rcpp::depends(RcppArmadillo)]]

//#include <gperftools/profiler.h>
#include <RcppArmadillo.h>
#include <vector>
#include <array>
#include <stdexcept>
#include <omp.h>
#include <unordered_map>
#include <iostream>
#include <string>
#include <cmath>
#include <boost/functional/hash.hpp>
#include <sys/resource.h>

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
    //mat * alignmentBinList ;
    std::vector<mat> alignmentBinList ;
    std::vector< std::vector<mat> > multiAlignmentBinList ;
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
    //phylo(const Rcpp::NumericMatrix &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::NumericVector, const bool, const int, const bool, const Rcpp::NumericVector &) ;
    phylo(const NumericMatrix &, const NumericVector &, const List &, const int, const bool, const bool, const uvec &) ;
    phylo() ;
    void logLikPhylo(const bool) ;
    std::vector<mat> getAlignmentBin(const uint) ;
    umat getAlignment() ;
    uint getNnode() ;
    uint getNumStates() ;
    uint getNumRateCats() ;
    uint getNumTips() ;
    SEXP getLogLik() ;
    SEXP getAlignmentBinSEXP() ;
    SEXP getPruningMat() ;
    //phylogeny() : numOpenMPthreads(4) {} ;
} ;

class phylogenyAlpha: public phylo {
public:
    // constructor, sets up data structures
    //phylogenyAlpha(const Rcpp::NumericMatrix &, const Rcpp::CharacterMatrix &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::CharacterVector &, const Rcpp::NumericVector, const bool, const bool, const int, const NumericVector &);
    phylogenyAlpha(const NumericMatrix &, const CharacterMatrix &, const NumericVector &, const List &, const int, const CharacterVector &, const bool, const bool, const uvec &) ;
    //phylogenyAlpha(const Rcpp::NumericMatrix &, const Rcpp::List &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::NumericVector, const bool, const bool, const int, const NumericVector &) ; // This is the faster constructor, I guess, since it doesn't need to perform alignment conversion.
    phylogenyAlpha(const NumericMatrix &, const List &, const NumericVector &, const List &, const int, const bool, const bool, const uvec &) ;
    phylogenyAlpha() ; // Empty constructor
private:
    Rcpp::CharacterMatrix alignmentRcpp ;
    void convertSTL() ;
    void convertSeqAlpha() ;
    void convertNucleoToNum(uint locusIndex) ;
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
    logTransMatVec.resize(numRateCats) ;
    sitePatterns = sitePatternsVec ;
    numLoci = sitePatterns.n_rows ;
    numUniqueLoci = sitePatterns.max() ;

    for (int i = 0; i < numRateCats ; i++) {

        logTransMatVec[i] = Rcpp::as<mat>(logTransMatList[i]) ;
    }

    children.resize(edge.n_rows + 2) ;

    for(uint i = 0; i < edge.n_rows; i++)
    {

        children[edge(i, 0)].push_back(edge(i, 1));
    }

    //childNodeInClusIndic = Rcpp::as<uvec>(RcppChildInClusIndic) ;
    if (returnMatIndic) {

        pruningMatVec.resize(numRateCats) ;
        for (uint i = 0 ; i < numRateCats; i++) {

            pruningMatVec[i] = Rcpp::NumericMatrix(numStates, numLoci) ;
        }
    } else{}
    if (internalPhyloFlag) {

        multiAlignmentBinList.resize(numRateCats) ;
        for (uint i = 0; i < numRateCats; i++) {

            multiAlignmentBinList[i].resize(numUniqueLoci) ;
        }
    }
}

phylogenyAlpha::phylogenyAlpha(const NumericMatrix & edgeMat, const CharacterMatrix & alignmentAlphaMat, const NumericVector & logLimProbsVec, const List & logTransMatList, const int numOpenMP, const CharacterVector & equivVec, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) : phylo(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsVec) {

    alignmentRcpp = alignmentAlphaMat ;
    equivalency = Rcpp::as<std::vector<std::string>>(equivVec) ;
    convertSTL() ;
    defineMap() ;
    convertSeqAlpha() ;
    if (!(logLimProbsVec[0] == 1000.0)) { // For getConvertedAlignmentToWrao, edgeMat is a placeholder for a phylo with two tips, but alignmentAlphaMat.nrow() is not 2. This is why this last part won't run. I use a placeholder for logLimProbsVec too, with first element equal to 1000, hence this check.
      numTips = alignmentAlphaMat.nrow() ;
      Nnode = edgeMat.nrow() - numTips + 1 ;
      compUpdateVec() ;
    } else{}
}

// This is the faster constructor, that does not convert the character alignment into matrices of identity vectors, because it's provided by the user (alignmentList).
phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::List & alignmentList, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numOpenMP, const bool returnMatIndic, const bool internalFlag, const uvec & sitePatternsVec) : phylo(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsVec) {

    if (!internalFlag) {

        numTips = Rcpp::as<mat>(alignmentList[0]).n_cols ;
        numUniqueLoci = alignmentList.size() ;
        alignmentBinList.resize(numUniqueLoci) ;
        for (int i = 0; i < numUniqueLoci ; i++) {

            alignmentBinList[i] = Rcpp::as<mat>(alignmentList[i]) ;
        }
    } else { // Tip vectors depend on the rate category.

        numTips = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(alignmentList[0])[0]).n_cols ;
        numUniqueLoci = Rcpp::as<Rcpp::List>(alignmentList[1]).size() ;
        for (uint i = 0; i < numRateCats ; i++) {

            for (int j = 0; j < numUniqueLoci ; j++) {

                multiAlignmentBinList[i][j] = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(alignmentList[i])[j]) ;
            }
        }
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
    for (uint i = 0; (i < Nnode) ; i++) {

        numChildren(i) = children[i+numTips+1].size() ;
    }
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

void phylogenyAlpha::convertSeqAlpha() {

    alignmentBinList.resize(numUniqueLoci) ;
    for (uint i = 0; i < numUniqueLoci; i++) {
        convertNucleoToNum(i) ;
    }
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

void phylogenyAlpha::convertNucleoToNum(uint locusIndex)  {

    mat containerMat(numStates, alignmentRcpp.nrow(), fill::zeros) ;
    std::vector<std::string> locusCol = alignmentAlpha[locusIndex] ;
    for (uint i = 0; i < containerMat.n_cols; i++) {

        containerMat.col(i) = mymap[locusCol[i]] ;
    }
    alignmentBinList[locusIndex] = containerMat ;
}

// Functions to modify!

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

        nodeTipMat.cols(0, numTips-1) = log(alignmentBinList[locusNum]) ;
    } else {

        nodeTipMat.cols(0, numTips-1) = multiAlignmentBinList[rateIndex][locusNum] ; // The alignment now depends on the rate category.
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
// End: Functions to modify
std::vector<mat> phylo::getAlignmentBin(const uint pos) {

    return (alignmentBinList) ;
}

SEXP phylo::getAlignmentBinSEXP() {

    std::vector<Rcpp::NumericVector> containerVec ;
    containerVec.resize(numUniqueLoci) ;
    for (uint i = 0 ; i< numUniqueLoci ; i++) {

        Rcpp::NumericVector myVector = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(alignmentBinList[i])) ;
        myVector.attr("dim") = Rcpp::Dimension(alignmentBinList[i].n_rows, alignmentBinList[i].n_cols);
        containerVec[i] = myVector ;
    }
    return Rcpp::wrap(containerVec) ;
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

NumericMatrix convMatToRcpp(arma::mat x) {
  NumericMatrix y = wrap(x) ;
  return(y) ;
}

SEXP phylo::getPruningMat() {

    std::vector<NumericMatrix> pruningMatVecReadj(pruningMatVec.size()) ;

    for (uint i = 0 ; i < pruningMatVec.size(); i++) {

      pruningMatVecReadj[i] = convMatToRcpp(as<mat>(pruningMatVec[i]).cols(sitePatterns - 1)) ;
    }
    return Rcpp::wrap(pruningMatVecReadj) ;
}

void phylo::logLikPhylo(const bool returnMatIndic) {

    mat logLikMat(numUniqueLoci, numRateCats, fill::ones);

    for (uint rateNum = 0; (rateNum < numRateCats); rateNum++) {

        #pragma omp parallel for
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
      numLoci = as<List>(alignmentBin[0]).size() ;
      sitePatternsAmended.resize(numLoci) ;
      sitePatternsAmended = conv_to<uvec>::from(linspace(1, numLoci, numLoci))  ;
    } else{

      sitePatternsAmended = as<uvec>(sitePatterns) ;
    }
    SEXP container ;
    //         ProfilerStart("/home/villandre/profileOut.out") ;

    phylogenyAlpha phyloObject(edgeMat, alignmentBin, logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsAmended) ;
    //#pragma omp parallel
    {
      phyloObject.logLikPhylo(returnMatIndic);
    }
    //             ProfilerStop();
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

    return phyloObject.getAlignmentBinSEXP() ;
}

// [[Rcpp::export]]

SEXP redimMultiBinByClus(Rcpp::List multiBinByClus) {

    uint numClus = multiBinByClus.size() ;
    uint numRateCats = (Rcpp::as<Rcpp::List>(multiBinByClus[0])).size() ;
    uint numStates = Rcpp::as<Rcpp::NumericMatrix>((Rcpp::as<Rcpp::List>(multiBinByClus[0]))[0]).nrow() ;
    uint numLoci = Rcpp::as<Rcpp::NumericMatrix>((Rcpp::as<Rcpp::List>(multiBinByClus[0]))[0]).ncol() ;
    std::vector<std::vector<Rcpp::NumericMatrix>> outputVec(numRateCats);

    for(uint i = 0; i < numRateCats ; i++) {

        outputVec.at(i).resize(numLoci) ;
        for (uint j = 0; j < numLoci; j++) {

            outputVec.at(i).at(j) = Rcpp::NumericMatrix(numStates, numClus) ;
            (outputVec.at(i).at(j)).attr( "dimnames" )  = Rcpp::List::create(R_NilValue, multiBinByClus.names()) ;
        }
    }

    for(uint i = 0; i < numLoci ; i++) {

        for (uint j = 0; j < numRateCats; j++) {

            for (uint k = 0; k < numClus; k++) {

               (outputVec.at(j).at(i))(Rcpp::_,k) = Rcpp::as<Rcpp::NumericMatrix>((Rcpp::as<Rcpp::List>(multiBinByClus[k]))[j])(Rcpp::_, i) ;
            }
        }
    }
    return(Rcpp::wrap(outputVec)) ;
}
// [[Rcpp::export]]

SEXP logLikCppToWrapV(List & edgeMatList, NumericVector & logLimProbsVec, List & logTransMatList, int numOpenMP, SEXP & equivVector, List alignmentBinList, const bool returnMatIndic, const bool internalFlag, const List sitePatternsList) {
//     #pragma omp parallel
//     {
  omp_set_num_threads(numOpenMP) ;

  uint numEvals = edgeMatList.size() ;
  IntegerVector sitePatterns(as<IntegerVector>(sitePatternsList[0]).size()) ;  // The number of loci does not change across clusters, hence the hard-coded 0.

  std::vector<SEXP> container(numEvals) ;
  std::vector<phylogenyAlpha> phylogenyAlphaContainer(numEvals) ;

  uvec sitePatternsAmended ;
  uint numLoci ;
  //#pragma omp parallel for if(numEvals > 7)
  for (uint i = 0; i < numEvals; i++) {

    if (!internalFlag) {

      numLoci = as<List>(sitePatternsList[i]).size() ;
      sitePatternsAmended.resize(numLoci) ;
      sitePatternsAmended = as<uvec>(sitePatternsList[i]) ;
    } else {

      numLoci = as<List>(alignmentBinList[i]).size() ;
      sitePatternsAmended = conv_to<uvec>::from(conv_to<uvec>::from(linspace(1, numLoci, numLoci))) ; // Site patterns do not apply for internal phylogenetic computations. A placeholder is therefore created here.
    }

    phylogenyAlphaContainer[i] = phylogenyAlpha(Rcpp::as<NumericMatrix>(edgeMatList[i]), Rcpp::as<List>(alignmentBinList[i]), logLimProbsVec, logTransMatList, numOpenMP, returnMatIndic, internalFlag, sitePatternsAmended) ;
  }

  for (int i = 0; i < numEvals; i++) {

    //ProfilerStart("/home/villandre/profileOut.out") ;
    phylogenyAlphaContainer[i].logLikPhylo(returnMatIndic) ;
    //ProfilerStop();
  }

  for (uint i = 0; i < numEvals; i++) {

    if (returnMatIndic) {

      container[i] = phylogenyAlphaContainer[i].getPruningMat() ;
    } else {

      container[i] = phylogenyAlphaContainer[i].getLogLik();
    }
  }
  return Rcpp::wrap(container) ;
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
