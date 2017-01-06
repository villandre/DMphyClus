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
    //umat alignment ;
    std::vector<mat> transMatVec ;
    vec limProbs ;
    uvec updateOrder ;
    double logLik ;
    std::vector<Rcpp::NumericMatrix> pruningMatVec ; // Each element of pruningMatVec comprises all the vectors at the root node used in computing the likelihood with the pruning algorithm. We get the likelihood for one rate category doing sum(log(pruningMat[[x]]%*%limProbs)). It follows that pruningMat[[x]] has one column per locus and one row per state (4 for nucleotide alignments, 20 for AA). This vector has as many elements as rate categories.
    //mat * alignmentBinList ;
    std::vector<mat> alignmentBinList ;
    std::vector< std::vector<mat> > multiAlignmentBinList ;
    uint numStates;
    uint numRateCats;
    std::vector<std::vector<int> > children;
    int numOpenMPthreads ;
    uvec childNodeInClusIndic ; // This is a vector with 0s and 1s. 1 means that the child node is at the end of an edge inside the cluster. 0 means that the child node is at the end of an edge outside a cluster. This is probably not functional and should be removed, until it's properly implemented.
    bool internalPhyloFlag ; // This flag indicates that this is an internal phylogeny, i.e. a phylogeny that supports subphylogenies for clusters. Concretely, this implies that the configuration corresponding to cluster centroids, that are at some of the tips of this phylogeny, are vectors with as many elements as rate categories.

    // Member functions:
    // Compute the order in which nodes must be update using Felsenstein's tree-pruning algorithm.
    void compUpdateVec() ;
    // Output node children:
    uvec Children(const umat &, const uint) ;
    // Member function: The likelihood at one locus for one rate category.
    double likOneLocusOneRate(const uint, const int, const bool) ;
    inline vec getLikVec(const uint, mat &, const int) ;
    inline void internalFun(const uint, mat &, const int) ;

public:
    // constructor, sets up data structures
    phylo(const Rcpp::NumericMatrix &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::NumericVector, const bool, const int, const bool) ;
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
    phylogenyAlpha(const Rcpp::NumericMatrix &, const Rcpp::CharacterMatrix &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::CharacterVector &, const Rcpp::NumericVector, const bool, const bool, const int);
    phylogenyAlpha(const Rcpp::NumericMatrix &, const Rcpp::List &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::NumericVector, const bool, const bool, const int) ; // This is the faster constructor, I guess, since it doesn't need to perform alignment conversion.
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

phylo::phylo(const Rcpp::NumericMatrix & edgeMat, const Rcpp::NumericVector & limProbsVec, const Rcpp::List & transMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::NumericVector RcppChildInClusIndic, const bool returnMatIndic, const int numOfLoci, const bool internalFlag) {

    mat edgeDouble = Rcpp::as<mat>(edgeMat); // The first argument (matrix) of the R function is called edgeMatrix.
    edge = conv_to<umat>::from(edgeDouble); // edge is a matrix of unsigned integers.
    limProbs = Rcpp::as<vec>(limProbsVec); // The second argument (numeric) of the R function is called limProbs.
    Nnode = NnodeCons ;
    numRateCats = numRateCatsCons ;
    numLoci = numOfLoci ;
    numStates = numStatesCons ;
    numOpenMPthreads = numOpenMP ;
    internalPhyloFlag = internalFlag ;
    numTips = edge.n_rows + 1 - Nnode ;
    transMatVec.resize(numRateCats) ;

    for (int i = 0; i < transMatList.size() ; i++) {

        transMatVec[i] = Rcpp::as<mat>(transMatList[i]) ;
    }
    children.resize(Nnode + numTips + 1) ;

    for(uint i = 0; i < edge.n_rows; i++)
    {
        children[edge(i, 0)].push_back(edge(i, 1));
    }
    compUpdateVec() ;
    childNodeInClusIndic = Rcpp::as<uvec>(RcppChildInClusIndic) ;
    if (returnMatIndic) {

        pruningMatVec.resize(numRateCats) ;
        for (uint i = 0 ; i < numRateCats; i++) {

            pruningMatVec[i] = Rcpp::NumericMatrix(numStates, numLoci) ;
        }
    } else{}
    if (internalPhyloFlag) {

        multiAlignmentBinList.resize(numRateCats) ;
        for (uint i = 0; i < numRateCats; i++) {

            multiAlignmentBinList[i].resize(numLoci) ;
        }
    }
}

phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::CharacterMatrix & alignmentAlphaMat, const Rcpp::NumericVector & limProbsVec, const Rcpp::List & transMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::CharacterVector & equivVec, const Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const int numLoci) : phylo(edgeMat, limProbsVec, transMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, numLoci, internalFlag) {

    alignmentRcpp = alignmentAlphaMat ;
    equivalency = Rcpp::as<std::vector<std::string>>(equivVec) ;
    convertSTL() ;
    defineMap() ;
    convertSeqAlpha() ;
    numTips = alignmentAlphaMat.nrow() ;
}

// This is the faster constructor, that does not convert the character alignment into matrices of identity vectors, because it's provided by the user (alignmentList).
phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::List & alignmentList, const Rcpp::NumericVector & limProbsVec, const Rcpp::List & transMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const int numLoci) : phylo(edgeMat, limProbsVec, transMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, numLoci, internalFlag) {

    if (!internalFlag) {

        alignmentBinList.resize(numLoci) ;
        for (int i = 0; i < numLoci ; i++) {

            alignmentBinList[i] = Rcpp::as<mat>(alignmentList[i]) ;
        }
    } else { // Tip vectors depend on the rate category.

        for (uint i = 0; i < numRateCats ; i++) {

            for (int j = 0; j < numLoci ; j++) {

                multiAlignmentBinList[i][j] = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(alignmentList[i])[j]) ;
            }
        }
    }
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

    uint numNodesTotal = numTips + Nnode ;
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

            updateOrderVec(indexForUpdateVec) = numTips + i + 1; // We want an updateOrderVec that reflects the one produced by R.
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

void phylogeny::convertSeq() {

    alignmentBinList.resize(numLoci) ;
    uvec firstOffset = cumsum(uvec(numTips).fill(numStates)) ;

    for (uint i = 0; (i < numLoci); i++) {

        mat binMatrix(numStates, numTips, fill::zeros) ; // This matrix has as many rows as sequences in the alignment, under the assumption that the alignment has one row per sequence.
        uvec oneLocusConf = alignment.col(i) ;
        uvec matIndices = oneLocusConf + firstOffset ;
        matIndices(span::all) += (-(numStates + 1)) ;
        binMatrix.elem(matIndices).ones() ;
        alignmentBinList[i] = binMatrix ;
    }
}

void phylogenyAlpha::convertSeqAlpha() {

    alignmentBinList.resize(numLoci) ;
    for (uint i = 0; i < numLoci; i++) {
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
inline vec phylo::getLikVec(const uint childNodeIndex, mat & nodeTipMat, const int rateIndex) {

    vec likVec((nodeTipMat).n_rows) ;
    likVec = (transMatVec[rateIndex])*((nodeTipMat).col(childNodeIndex - 1)) ; // Indices start at 0, not at 1.
    return likVec ;
}

inline void phylo::internalFun(const uint parentNum, mat & nodeTipMat, const int rateIndex) {

    vec likVecByChild(nodeTipMat.n_rows, fill::ones) ;
    uint numChildren = children[parentNum].size() ;
    uvec theChildren(numChildren) ;
    for(uint i = 0; i < numChildren; i++) {

        theChildren(i) = children[parentNum][i];
    }
    for (uint i = 0 ; (i < numChildren); i++) {

        likVecByChild = likVecByChild % getLikVec(theChildren(i), nodeTipMat, rateIndex) ;
    }

    nodeTipMat.col(parentNum - 1) = likVecByChild ;
}

double phylo::likOneLocusOneRate(const uint locusNum, const int rateIndex, const bool returnMatIndic) {

    mat nodeTipMat = zeros<mat>(numStates, edge.max()) ; // This is a matrix that stores likelihood vectors. Each vector has as many elements as potential states.

    if (!internalPhyloFlag) {

        nodeTipMat.cols(0, numTips-1) = alignmentBinList[locusNum] ;
        nodeTipMat = nodeTipMat*pow(10, (double) 300/ (double) numTips) ;// To prevent the case where values fall below the minimum allowed by doubles. This scaling is taken into account at the end, when a log-lik value is computed.
    } else {

        nodeTipMat.cols(0, numTips-1) = multiAlignmentBinList[rateIndex][locusNum] ; // The alignment now depends on the rate category.
        double myFactor =  (double) 300 / double (numTips) - (double) 300 ;
        nodeTipMat = nodeTipMat * pow(10, myFactor) ;// To prevent the case where values fall below the minimum allowed by doubles. This scaling is taken into account at the end, when a log-lik value is computed.
    }

    for (uint i = 0; (i < updateOrder.n_rows) ; i++) {

        internalFun(updateOrder(i), nodeTipMat, rateIndex) ;
    }
    double likForOneLocusOneRate = sum(nodeTipMat.col(numTips)%limProbs) ; // No +1 since I start at 0 instead of 1.

    if (returnMatIndic) {
        for (uint i = 0 ; i < nodeTipMat.n_rows; i++) {
            double valueInMatrix = nodeTipMat.at(i,numTips) ;
            pruningMatVec[rateIndex](i, locusNum) = valueInMatrix ;
        }
    } else{}

    return likForOneLocusOneRate ;
}
// End: Functions to modify
std::vector<mat> phylo::getAlignmentBin(const uint pos) {

    return (alignmentBinList) ;
}

SEXP phylo::getAlignmentBinSEXP() {

    std::vector<Rcpp::NumericVector> containerVec ;
    containerVec.resize(numLoci) ;
    for (uint i = 0 ; i< numLoci ; i++) {

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

SEXP phylo::getPruningMat() {

//     std::vector<Rcpp::NumericMatrix> pruningMatFixed(pruningMatVec.size()) ;
//     for (uint i = 0 ; i < pruningMatFixed.size() ; i++) {
//         pruningMatFixed.at(i) = pruningMatVec.at(i)*pow(10, -300) ;
//     }
//     return Rcpp::wrap(pruningMatFixed) ;
    return Rcpp::wrap(pruningMatVec) ; // We can ignore the scaling as long as there's an offset when we compute the log-lik.
}

void phylo::logLikPhylo(const bool returnMatIndic) {

    mat likMat(numLoci, numRateCats, fill::ones);
//    #pragma omp parallel for
    for (uint rateNum = 0; (rateNum < numRateCats); rateNum++) {

         #pragma omp parallel for
        for(uint locusNum = 0; locusNum < numLoci; locusNum++) {

            likMat.at(locusNum, rateNum) = likOneLocusOneRate(locusNum, rateNum, returnMatIndic) ; // We applied a scaling factor that will have to be fixed.
        }
    }
    vec logLikSumRate = log(sum(likMat, 1)) - (double) 300 * log(10); // The offset compensates for the scaling.
    vec meanEffectVec = vec(numLoci).fill(log(numRateCats)) ;
    logLik = sum(logLikSumRate - meanEffectVec)  ;
}

template<class Mat>
void print_matrix(Mat matrix) {

    //matrix.print(std::cout);
}


template<class Col>
void print_vector(Col colvec) {

    //colvec.print(std::cout);
}

//provide explicit instantiations of the template function for
//every matrix type you use somewhere in your program.
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);
template void print_vector<arma::uvec>(arma::uvec colvec);
template void print_vector<arma::vec>(arma::vec colvec);

// [[Rcpp::export]]

SEXP logLikCppToWrap(NumericMatrix & edgeMat, NumericVector & limProbsVec, List & transMatList, int numOpenMP, SEXP & equivVector, List alignmentBin, NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const NumericVector sitePatterns) {

    omp_set_num_threads(numOpenMP) ;
    uint numTips ;
    uint numLoci ;
    SEXP container ;
    uint numStatesCons = limProbsVec.size() ;
    uint numRateCatsCons = transMatList.size() ;
    uint NnodeCons ;

    if (internalFlag) {

        numTips = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(alignmentBin[0])[0]).n_cols ; 
        numLoci = (Rcpp::as<Rcpp::List>(alignmentBin[0])).size() ;
    } else {

        numTips = Rcpp::as<mat>(alignmentBin[0]).n_cols ;
        numLoci = alignmentBin.size() ;
    }
    
    NnodeCons = edgeMat.nrow() - numTips + 1 ;
    //         ProfilerStart("/home/villandre/profileOut.out") ;
    
    
    phylogenyAlpha phyloObject(edgeMat, alignmentBin, limProbsVec, transMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, internalFlag, numLoci) ;
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

SEXP getConvertedAlignmentToWrap(uint numStatesCons, int numOpenMP, SEXP & equivVector, Rcpp::CharacterMatrix & alignmentAlphaMat) {

    double src[] = {3, 3, 1, 2};
    Rcpp::NumericMatrix placeholderMat(2, 2, src) ; //src is considered an iterator.
    Rcpp::NumericVector placeholderVec = Rcpp::NumericVector::create(1000,1000) ;
    Rcpp::List placeholderList(1) ; // This is a placeholder: 1 is an arbitrary number.
    placeholderList[0] = Rcpp::NumericMatrix(2,2) ;
    bool returnMatIndic = false ;
    bool internalFlag = false ;
    int placeholderNnodeCons = 1;
    int placeholderNumRateCats = 4;

    uint numLoci = alignmentAlphaMat.ncol() ;
    phylogenyAlpha phyloObject = phylogenyAlpha(placeholderMat, alignmentAlphaMat, placeholderVec, placeholderList, numStatesCons, placeholderNumRateCats, placeholderNnodeCons, numOpenMP, equivVector, placeholderVec, returnMatIndic, internalFlag, numLoci) ;

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

SEXP logLikCppToWrapV(List & edgeMatList, NumericVector & limProbsVec, List & transMatList, int numOpenMP, SEXP & equivVector, List alignmentBinList, NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const List sitePatternsList) {
//     #pragma omp parallel
//     {
  omp_set_num_threads(numOpenMP) ;
  
  uint numEvals = edgeMatList.size() ;
  uint numStatesCons = limProbsVec.size() ;
  uint numRateCatsCons, numTips, numBranches, NnodeCons, numLoci ;

  std::vector<SEXP> container(numEvals) ;
  std::vector<phylogenyAlpha> phylogenyAlphaContainer(numEvals) ;
  
  Rcpp::List transMatListOneSize ;
  
  numRateCatsCons = transMatList.size() ;
  
  if (internalFlag) {
    
    numLoci = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(alignmentBinList[0])[0]).size() ;
  } else {
    
    numLoci = Rcpp::as<Rcpp::List>(alignmentBinList[0]).size() ;
  }
  
  //#pragma omp parallel for if(numEvals > 7)
  for (uint i = 0; i < numEvals; i++) {
    numBranches  = Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i]).nrow() ;
    numTips = Rcpp::min(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i])(Rcpp::_,0)) - 1 ; 
    NnodeCons = Rcpp::max(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i])(Rcpp::_,0)) - numTips ; 
    
    transMatListOneSize = transMatList ;
    
    phylogenyAlphaContainer[i] = phylogenyAlpha(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i]), Rcpp::as<Rcpp::List>(alignmentBinList[i]), limProbsVec, transMatListOneSize, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, internalFlag, numLoci) ;
  }
  
  for (int i = 0; i < edgeMatList.size(); i++) {
    
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

// template <class Container>
// std::vector<std::size_t> get_hashes(Container const& x)
// {
//   std::vector<std::size_t> hashes;
//   std::transform(x.begin(), x.end(), hashes.begin(), boost::hash<typename Container::value_type>());
//   
//   return hashes;
// }
// 
// // [[Rcpp::export]]
// 
// Rcpp::NumericVector getSitePatterns(Rcpp::CharacterMatrix alignmentMat) {
// 
//   std::vector<Col> sitePatternLetters ;
//   std::vector<std::vector<string>> sitePatternHash ;
//   sitePatterns = get_hashes(alignmentMat) ;
//   return Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(sitePatterns)) ;
// }



