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

//#include <boost/algorithm/string.hpp>
// [[Rcpp::plugins(openmp)]]

using namespace arma;

class phylo { // NOTE: ONCE CODE IS DEBUGGED, YOU CAN SPEED IT UP BY USING [] or .at() FOR INDEXING.
protected:
    umat edge ;
    uint Nnode ;
    uint numTips ;
    uint numLoci ;
    //umat alignment ;
    std::vector<mat> logTransMatVec ;
    vec logLimProbs ;
    uvec updateOrder ;
    double logLik ;
    std::vector<Rcpp::NumericMatrix> pruningMatVec ; // Each element of pruningMatVec comprises all the vectors at the root node used in computing the likelihood with the pruning algorithm. We get the likelihood for one rate category doing sum(log(pruningMat[[x]]%*%logLimProbs)). It follows that pruningMat[[x]] has one column per locus and one row per state (4 for nucleotide alignments, 20 for AA). This vector has as many elements as rate categories.
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
    double logLikOneLocusOneRate(const uint, const int, const bool) ;
    vec getLogLikVec(const uint, mat &, const int) ;
    void internalFun(const uint, mat &, const int) ;

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

class phylogeny: public phylo {
public:
    // constructor, sets up data structures
    phylogeny(const Rcpp::NumericMatrix &, const Rcpp::NumericMatrix &, const Rcpp::NumericVector &, const Rcpp::List &, const int, const int, const int, const int, const Rcpp::NumericVector, const bool) ;
    phylogeny() ; // Empty constructor
private:
    void convertSeq() ;
    umat alignment ;
 };

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

phylo::phylo(const Rcpp::NumericMatrix & edgeMat, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::NumericVector RcppChildInClusIndic, const bool returnMatIndic, const int numOfLoci, const bool internalFlag) {

    mat edgeDouble = Rcpp::as<mat>(edgeMat); // The first argument (matrix) of the R function is called edgeMatrix.
    edge = conv_to<umat>::from(edgeDouble); // edge is a matrix of unsigned integers.
    logLimProbs = Rcpp::as<vec>(logLimProbsVec); // The second argument (numeric) of the R function is called logLimProbs.
    Nnode = NnodeCons ;
    numRateCats = numRateCatsCons ;
    numLoci = numOfLoci ;
    numStates = numStatesCons ;
    numOpenMPthreads = numOpenMP ;
    internalPhyloFlag = internalFlag ;
    numTips = edge.n_rows + 1 - Nnode ;
    logTransMatVec.resize(numRateCats) ;

    for (int i = 0; i < logTransMatList.size() ; i++) {

        logTransMatVec[i] = Rcpp::as<mat>(logTransMatList[i]) ;
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

phylogeny::phylogeny(const Rcpp::NumericMatrix & edgeMat, const Rcpp::NumericMatrix & alignmentMat, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic) : phylo(edgeMat, logLimProbsVec, logTransMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, alignmentMat.ncol(), returnMatIndic) {

    mat alignmentDouble = Rcpp::as<mat>(alignmentMat); // The third argument (matrix) of the R function is called alignment, and so on...
    numLoci = alignmentDouble.n_cols ;
    alignment = conv_to<umat>::from(alignmentDouble); // edge is a matrix of unsigned integers.
    numTips = alignment.n_rows ;
    convertSeq() ;
}

phylogeny::phylogeny() : phylo() {} //Mostly for setting up container vectors. Doesn't seem very efficient though!

phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::CharacterMatrix & alignmentAlphaMat, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::CharacterVector & equivVec, const Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const int numLoci) : phylo(edgeMat, logLimProbsVec, logTransMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, numLoci, internalFlag) {

    alignmentRcpp = alignmentAlphaMat ;
    equivalency = Rcpp::as<std::vector<std::string>>(equivVec) ;
    convertSTL() ;
    defineMap() ;
    convertSeqAlpha() ;
    numTips = alignmentAlphaMat.nrow() ;
}

// This is the faster constructor, that does not convert the character alignment into matrices of identity vectors, because it's provided by the user (alignmentList).
phylogenyAlpha::phylogenyAlpha(const Rcpp::NumericMatrix & edgeMat, const Rcpp::List & alignmentList, const Rcpp::NumericVector & logLimProbsVec, const Rcpp::List & logTransMatList, const int numStatesCons, const int numRateCatsCons, const int NnodeCons, const int numOpenMP, const Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const int numLoci) : phylo(edgeMat, logLimProbsVec, logTransMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, numLoci, internalFlag) {

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
    
    vec logLikVecByChild(nodeTipMat.n_rows, fill::zeros) ;
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

    mat logLikMat(numLoci, numRateCats, fill::ones);
//    #pragma omp parallel for
    for (uint rateNum = 0; (rateNum < numRateCats); rateNum++) {
        
         #pragma omp parallel for
        for(uint locusNum = 0; locusNum < numLoci; locusNum++) {

            logLikMat.at(locusNum, rateNum) = logLikOneLocusOneRate(locusNum, rateNum, returnMatIndic) ; 
        }
    }
    vec logLikSumRate = sum(logLikMat, 1) ;
    vec meanEffectVec = vec(numLoci).fill(log(numRateCats)) ;
    logLik = sum(logLikSumRate - meanEffectVec)  ;
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

SEXP logLikCppToWrap(Rcpp::NumericMatrix & edgeMat, Rcpp::NumericMatrix & alignmentMat, Rcpp::NumericVector & logLimProbsVec, Rcpp::List & logTransMatList, int numOpenMP, SEXP & equivVector, Rcpp::CharacterMatrix & alignmentAlphaMat, Rcpp::List alignmentBin, Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag) {
    
    //omp_set_dynamic(0) ; // Is it required? No... It seems to be a default setting.
    omp_set_num_threads(numOpenMP) ;
    bool testBool ;
    uint numLoci ;
    //uint numTips ;
    mat firstMatrix ; // Is firstMatrix useful?
    SEXP container ;
    uint numStatesCons = logLimProbsVec.size() ;
    uint numRateCatsCons = logTransMatList.size() ;
    uint NnodeCons ;

    if (internalFlag) {

        firstMatrix = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(alignmentBin[0])[0]) ;
        numLoci = (Rcpp::as<Rcpp::List>(alignmentBin[0])).size() ;
    } else {

        firstMatrix = Rcpp::as<mat>(alignmentBin[0]) ;
        numLoci = alignmentBin.size() ;
    }
    testBool = firstMatrix[0,0] == 1000 ;
    //numTips = firstMatrix.n_cols ; // alignmentBin may not be defined, in which case a placeholder is fed to the function. This is why , if testBool is TRUE, we also need to update numTips.

    if (testBool) { // 1000 is a placeholder indicating missingness. Is there a better way to do this?

        if (alignmentMat(0,0) == 0) { // Again, in that case, we have a placeholder for alignmentMat, which means that alignmentAlphaMat should be defined.
            //numTips = alignmentAlphaMat.nrow() ;
            NnodeCons = edgeMat.nrow() - alignmentAlphaMat.nrow() + 1; // Make sure this works.

            phylogenyAlpha phyloObject(edgeMat, alignmentAlphaMat, logLimProbsVec, logTransMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, equivVector, childNodeInClusIndic, returnMatIndic, internalFlag, numLoci) ;

            phyloObject.logLikPhylo(returnMatIndic);

            if (returnMatIndic) {

                container = phyloObject.getPruningMat();
            } else {

                container = phyloObject.getLogLik();
            }
        } else {

            NnodeCons = edgeMat.nrow() - alignmentMat.nrow() + 1 ;

            phylogeny phyloObject = phylogeny(edgeMat, alignmentMat, logLimProbsVec, logTransMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic) ;
            phyloObject.logLikPhylo(returnMatIndic) ;
            if (returnMatIndic) {

                container = phyloObject.getPruningMat() ;
            } else {

                container = phyloObject.getLogLik() ;
            }
        }
    } else {
//         ProfilerStart("/home/villandre/profileOut.out") ;
        NnodeCons = edgeMat.nrow() - firstMatrix.n_cols + 1; // Make sure this works.
        
        phylogenyAlpha phyloObject(edgeMat, alignmentBin, logLimProbsVec, logTransMatList, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, internalFlag, numLoci) ;
        //#pragma omp parallel
        
        phyloObject.logLikPhylo(returnMatIndic);
        
//             ProfilerStop();
        if (returnMatIndic) {

            container = phyloObject.getPruningMat();
        } else {

            container = phyloObject.getLogLik();
        }
    }

    return container ;
}

// [[Rcpp::export]]

SEXP getConvertedAlignmentToWrap(Rcpp::NumericMatrix & alignmentMat, uint numStatesCons, int numOpenMP, SEXP & equivVector, Rcpp::CharacterMatrix & alignmentAlphaMat) {

    double src[] = {3, 3, 1, 2};
    Rcpp::NumericMatrix placeholderMat(2, 2, src) ; //src is considered an iterator.
    Rcpp::NumericVector placeholderVec = Rcpp::NumericVector::create(1000,1000) ;
    Rcpp::List placeholderList(1) ; // This is a placeholder: 1 is an arbitrary number.
    placeholderList[0] = Rcpp::NumericMatrix(2,2) ;
    bool returnMatIndic = false ;
    bool internalFlag = false ;
    int placeholderNnodeCons = 1;
    int placeholderNumRateCats = 4;

    if (alignmentMat(0,0) == 0) {

        uint numLoci = alignmentAlphaMat.ncol() ;
        phylogenyAlpha phyloObject = phylogenyAlpha(placeholderMat, alignmentAlphaMat, placeholderVec, placeholderList, numStatesCons, placeholderNumRateCats, placeholderNnodeCons, numOpenMP, equivVector, placeholderVec, returnMatIndic, internalFlag, numLoci) ;

        return phyloObject.getAlignmentBinSEXP() ;
    } else {

        phylogeny phyloObject = phylogeny(placeholderMat, alignmentMat, placeholderVec, placeholderList, numStatesCons, placeholderNumRateCats, placeholderNnodeCons, numOpenMP, placeholderVec, returnMatIndic) ;
        return phyloObject.getAlignmentBinSEXP();
    }
}

// When we consider tips resulting from pruning a subtree, their configurations become probabilistic and so, depends on the number of rate variation categories. We must therefore create one alignmentBin per rate category.
//pruningMatList has as many elements as seqNumVec, that is, one element per tip to substitute. Each element of pruningMatList is itself a list with as many elements as rate categories. Finally, each element of the sub-list is a matrix with as many rows as states, and as many columns as loci.

// [[Rcpp::export]]aaaa

// SEXP amendBinAlignToWrap(const Rcpp::List & alignmentBinList, const Rcpp::NumericVector seqNumVec, const Rcpp::List pruningMatList, int numOpenMPthreads) {
//
//     //First, we recast the list of pruning matrices as a std::vector of vectors of matrices. The operation is parallelized. We then do the same for alignmentBinList.
//     uint numRateCats = pruningMatList.size() ;
//     uint numTipsToChange = seqNumVec.size() ;
//     uint numLoci = alignmentBinList.size() ;
//     uint lengthSubPruningList = Rcpp::as<Rcpp::List>(pruningMatList[0]).size() ;
//     if (numTipsToChange != lengthSubPruningList) {
//         throw std::invalid_argument("The number of elements in seqNumVec should match that in pruningMatList.\n");
//     }  else{}
//     std::vector< std::vector<Rcpp::NumericMatrix> > pruningMatVec(numRateCats) ; // Outer vector should have as many elements as rate categories. Inner vector should have as many elements as the number of tips to change.
//     std::vector< std::vector<Rcpp::NumericMatrix> > alignmentBinVec(numRateCats) ; // Outer vector should have as many elements as rate categories. Inner vector should have as many elements as loci.
//
//     for (uint i = 0; i < numRateCats; i++) {
//
//         pruningMatVec[i].resize(numTipsToChange) ;
//         for (uint j = 0; j < numTipsToChange; j++) {
//
//             Rcpp::NumericMatrix castedMatrix = Rcpp::as<Rcpp::List>(pruningMatList[i])[j] ;
//             pruningMatVec[i][j] = castedMatrix ;
//         }
//     }
//     std::vector<Rcpp::NumericMatrix> alignmentBinRecast(numLoci) ;
//
//     // Recasting alignmentBinList
//
//     for (uint i = 0; i < numLoci; i++) {
//
//         alignmentBinRecast.at(i) = Rcpp::as<Rcpp::NumericMatrix>(alignmentBinList[i]) ;
//     }
//
//     for (uint i = 0; i < numRateCats ; i++) {
//
//         alignmentBinVec[i].resize(numLoci) ;
//         alignmentBinVec[i] = alignmentBinRecast ;
//     } // This initializes alignmentBinVec. Each element has the same value and will need to be updated iteratively.
//
//     for (uint i = 0; i < numRateCats; i++) {
//
//         for (uint j = 0; j < numTipsToChange; j++) {
//
//             for (uint k = 0; k < numLoci; k++) {
//
//                 (alignmentBinVec[i][k])(Rcpp::_, seqNumVec[j] - 1) = (pruningMatVec[i][j])(Rcpp::_, k) ; // The -1 is there to account for R starting to count indices at 1 instead of 0.
//             }
//         }
//     }
//     return Rcpp::wrap(alignmentBinVec) ;
// }

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

SEXP logLikCppToWrapV(Rcpp::List & edgeMatList, Rcpp::List & alignmentMatList, Rcpp::NumericVector & logLimProbsVec, Rcpp::List & logTransMatList, int numOpenMP, SEXP & equivVector, Rcpp::List & alignmentAlphaMatList, Rcpp::List alignmentBinList, Rcpp::NumericVector childNodeInClusIndic, const bool returnMatIndic, const bool internalFlag, const bool priorBySizeTransMatBool) { // Function probably won't work with alignmentAlphaMatList...
//     #pragma omp parallel
//     {
    omp_set_num_threads(numOpenMP) ;
    bool alignmentBinListDefined ; // Indicates whether alignmentBinList is defined.
    uint numLoci ;

    mat firstMatrix ;
    uint numEvals = edgeMatList.size() ;
    std::vector<SEXP> container(numEvals) ;
    std::vector<phylogenyAlpha> phylogenyAlphaContainer(numEvals) ;
    std::vector<phylogeny> phylogenyContainer(numEvals) ;
    bool alphaMatPlaceholderFlag = Rcpp::as<Rcpp::CharacterMatrix>(alignmentAlphaMatList[0]).nrow() == 1;
    uint numStatesCons = logLimProbsVec.size() ;
    uint numRateCatsCons ;
    uint numTips ;
    uint numBranches ;
    uint NnodeCons ;

    if (priorBySizeTransMatBool) {

        numRateCatsCons = Rcpp::as<Rcpp::List>(logTransMatList[0]).size() ;
    } else {

        numRateCatsCons = logTransMatList.size() ;
    }
    Rcpp::List logTransMatListOneSize ;

    if (internalFlag) {

        firstMatrix = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(alignmentBinList[0])[0])[0]) ;
        numLoci = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(alignmentBinList[0])[0]).size() ;
        alignmentBinListDefined = true ; // For sure, in that case, alignmentBinList will be defined, and alignmentMatList will be a placeholder.

    } else {

        firstMatrix = Rcpp::as<mat>(Rcpp::as<Rcpp::List>(alignmentBinList[0])[0]) ;
        alignmentBinListDefined = firstMatrix[0,0] != 1000 ; // Determines if we have a placeholder for alignmentBinList (see the definition for logLikCppV in the R file).
        if (alignmentBinListDefined) {

            numLoci = Rcpp::as<Rcpp::List>(alignmentBinList[0]).size() ;
        } else {

            if (alphaMatPlaceholderFlag) {

                numLoci = Rcpp::as<Rcpp::NumericMatrix>(alignmentMatList[0]).ncol() ;
            } else {

                numLoci = Rcpp::as<Rcpp::CharacterMatrix>(alignmentAlphaMatList[0]).ncol() ;
            }
        }
    }

    // Create objects
    //#pragma omp parallel for if(numEvals > 7)
    for (uint i = 0; i < numEvals; i++) {
        numBranches  = Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i]).nrow() ;
        //numTips = numBranches/2 + 1 ; // This only works under the assumption that trees are bifurcating! This will have to be changed if the hypothesis is dropped.
        numTips = Rcpp::min(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i])(Rcpp::_,0)) - 1 ; // This works even when the tree is not bifurcating.
        //NnodeCons = numTips - 1 ; // This only works under the assumption that trees are bifurcating! This will have to be changed if the hypothesis is dropped.

        NnodeCons = Rcpp::max(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i])(Rcpp::_,0)) - numTips ; // This also works with multifurcating trees

        if (alignmentBinListDefined) { // 1000 is a placeholder indicating missingness. Is there a better way to do this?

            if (priorBySizeTransMatBool) {

                long long int numTipsRecast = numTips ;
                Rcpp::CharacterVector listNames(logTransMatList.size()) ;
                listNames = logTransMatList.names() ;
                std::string maxName = Rcpp::as<std::string>(listNames[logTransMatList.size()-1]) ;

                int maxSize = std::stoi(maxName) ;
                bool testSizeBool = (numTipsRecast <= maxSize) ;
                numTipsRecast = numTipsRecast*testSizeBool + maxSize*(!testSizeBool) ; // This assumes that the list elements are named and ordered, i.e. the last list element is for the largest cluster size in the distribution!

                logTransMatListOneSize = logTransMatList[std::to_string(numTipsRecast)] ;
            } else{

                logTransMatListOneSize = logTransMatList ;
            }

            phylogenyAlphaContainer[i] = phylogenyAlpha(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i]), Rcpp::as<Rcpp::List>(alignmentBinList[i]), logLimProbsVec, logTransMatListOneSize, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic, internalFlag, numLoci) ;
        } else {

            if (alphaMatPlaceholderFlag) { // We have a placeholder for alignmentAlphaMat.

                if (priorBySizeTransMatBool) {

                    long long int numTipsRecast = numTips ;
                    logTransMatListOneSize = logTransMatList[std::to_string(numTipsRecast)] ;
                } else{

                    logTransMatListOneSize = logTransMatList ;
                }

                phylogenyContainer[i] = phylogeny(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i]), Rcpp::as<Rcpp::NumericMatrix>(alignmentMatList[i]), logLimProbsVec, logTransMatListOneSize, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, childNodeInClusIndic, returnMatIndic) ;
            } else {

                if (priorBySizeTransMatBool) {

                    long long int numTipsRecast = numTips ;
                    logTransMatListOneSize = Rcpp::as<Rcpp::List>(logTransMatList[std::to_string(numTipsRecast)]) ; //Make sure the list labels carry over to c++
                } else{

                    logTransMatListOneSize = logTransMatList ;
                }

                phylogenyAlphaContainer[i] = phylogenyAlpha(Rcpp::as<Rcpp::NumericMatrix>(edgeMatList[i]), Rcpp::as<Rcpp::CharacterMatrix>(alignmentAlphaMatList[i]), logLimProbsVec, logTransMatListOneSize, numStatesCons, numRateCatsCons, NnodeCons, numOpenMP, equivVector, childNodeInClusIndic, returnMatIndic, internalFlag, numLoci) ;
            }
        }
    }

    for (int i = 0; i < edgeMatList.size(); i++) {

        if (alignmentBinListDefined) { // 1000 is a placeholder indicating missingness. Is there a better way to do this?

            //ProfilerStart("/home/villandre/profileOut.out") ;
            phylogenyAlphaContainer[i].logLikPhylo(returnMatIndic) ;
            //ProfilerStop();
        } else {

            if (alphaMatPlaceholderFlag) { // We have a placeholder for alignmentAlphaMatList

                phylogenyContainer[i].logLikPhylo(returnMatIndic) ;
            } else {

                phylogenyAlphaContainer[i].logLikPhylo(returnMatIndic) ;
            }
        }
    }
//     }

    for (uint i = 0; i < numEvals; i++) {

        if (alignmentBinListDefined) {

            if (returnMatIndic) {

                container[i] = phylogenyAlphaContainer[i].getPruningMat() ;
            } else {

                container[i] = phylogenyAlphaContainer[i].getLogLik();
            }

        } else {

            if (alphaMatPlaceholderFlag) { // Again, in that case, we have a placeholder for alignmentAlphaMat.

                if (returnMatIndic) {

                    container[i] = phylogenyContainer[i].getPruningMat();
                } else {

                    container[i] = phylogenyContainer[i].getLogLik();
                }
            } else {

               if (returnMatIndic) {

                    container[i] = phylogenyAlphaContainer[i].getPruningMat();
                } else {

                    container[i] = phylogenyAlphaContainer[i].getLogLik();
                }
            }
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


