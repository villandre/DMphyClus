#' Cluster an alignment with DM-PhyClus
#'
#' Implementation of the DM-PhyClus method described in Villandr\'{e} et al. 2016.
#'
#' @param numIters number of iterations for the MCMC sampler
#' @param numLikThreads number of openMP threads for log-likelihood evaluations
#' @param numMovesNNIbetween number of nearest-neighbour interchange moves
#' used in proposing transitions in the space of between-cluster phylogenies
#' @param numMovesNNIwithin number of nearest-neighbour interchange moves used
#' in proposing transitions in the space of within-cluster phylogenies
#' @param alignment matrix of characters. When using DNA data, all nucleotides
#' and ambiguities should be coded with small letters, following the IUPAC standard
#' @param startingValues list of starting values for the model parameters: elements
#'named phylogeny, clusInd, and alpha are required
#' @param numSamplesForTransMat number of random values used in Monte Carlo estimation
#' of the transition probability matrices
#' @param meanWithinBranchVec vector of positive numerical values giving a range of
#' potential mean branch lengths in the within-cluster phylogenies
#' @param meanBetweenBranchVec vector of positive numerical values giving a range of
#' potential mean branch lengths in the between-cluster phylogeny
#' @param limProbs vector of numerical values, with named elements, giving the limiting
#' probabilities for all states
#' @param clusPhyloUpdateProp number between 0 and 1 indicating the proportion of
#' within-cluster phylogenies to update in each MCMC iteration
#' @param numSplitMergeMoves the number of times the algorithm should try to split or
#' merge clusters in each MCMC iteration
#' @param numGammaCat number of among-loci substitution rate variation categories
#' @param discGammaPar value of the the discrete gamma distribution parameter tuning
#' among-sites rate variation
#' @param Qmatrix substitution rate matrix with rows and columns names corresponding
#' to states in the alignment
#' @param shapeForAlpha shape parameter of the gamma distribution used as a prior for
#' the Dirichlet concentration parameter
#' @param scaleForAlpha scale parameter of the gamma distribution used as a prior for
#' the Dirichlet concentration parameter
#' @param shiftForAlpha minimum value the Dirichlet concentration parameter can take
#' @param poisRateNumClus Poisson parameter for the number of clusters distribution;
#' if left unspecified, it is set equal to the number of clusters in the startingValues
#' @param betweenClusTransMatList list of lists of potential transition rate matrices
#' along branches in the supporting phylogeny. See details
#' @param withinClusTransMatList Same as betweenClusTransMatList, but for the
#' within-cluster phylogenies
#'
#' @details It is possible to supply transition probability matrices directly to the function
#' by specifying values for betweenClusTransMatList and withinClusTransMatList. Since estimating
#' those matrices might take a long time, it is preferable to perform the computations only once
#' using the \code{\link{outputTransMatList}} function, save the result to a file and re-use the
#' matrices whenever required. When those matrices are provided, no value need be specified for
#' Qmatrix, numGammaCat, discGammaPar, meanBetweenBranchVec, and meanWithinBranchVec.
#'
#' @return A list with two components:
#' \itemize{
#'  \item{chain:} {list with each element itself a list giving the sampled parameter values.
#' betweenTransMatListIndex and withinTransMatListIndex correspond to the index of the assumed mean
#' branch length in meanBetweenBranchVec and meanWithinBranchVec, respectively}
#'  \item{MAPestimate:}{ vector giving the maximum posterior probability cluster membership indices estimate}
#' }
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

DMphyClusChain <- function(numIters, numLikThreads = 1, numMovesNNIbetween = 1, numMovesNNIwithin = 1, alignment, startingValues, numSamplesForTransMat = 1e5, meanBetweenBranchVec, meanWithinBranchVec, limProbs, clusPhyloUpdateProp = 1, numSplitMergeMoves = 1, numGammaCat = 3, discGammaPar, Qmatrix, shapeForAlpha, scaleForAlpha, shiftForAlpha = 0, poisRateNumClus, betweenClusTransMatList, withinClusTransMatList) {

    .checkInput(startingValues = startingValues, Qmatrix = Qmatrix, alignment = alignment, limProbs = limProbs, shiftForAlpha = shiftForAlpha)
    if (!is.null(rownames(alignment))) { ## The tip ordering in the starting phylogeny should match the order of the rows in the alignment.
        alignment <- alignment[startingValues$phylogeny$tip.label,]
    } else{
        warning("The alignment has unnamed rows (sequences). Make sure that the order of the tips in startingValues$phylogeny matches that of the rows in alignment. \n")
        rownames(alignment) <- startingValues$phylogeny$tip.label
    }

    if (sum(limProbs) != 1) {
        limProbs <- limProbs/sum(limProbs)
        warning("limProbs does not sum to 1. Standardizing... \n")
    } else{}
    startingValues$clusInd <- startingValues$clusInd[startingValues$phylogeny$tip.label] ## Reordering the starting value for clusInd such that it matches the tips in the phylogeny...

    ## First, we obtain the transition rate matrices...
    if (missing(betweenClusTransMatList)) {
        QmatTest <- !is.null(rownames(Qmatrix))
        if (QmatTest) { ## We re-order the rows and columns of the Qmatrix, just in case...
            Qmatrix <- Qmatrix[names(limProbs),]
            Qmatrix <- Qmatrix[,names(limProbs)]
        } else{}
        cat("Estimating transition probability matrices for branches in the supporting phylogeny... ")
        allIntMatList <- lapply(meanBetweenBranchVec, FUN = function(x) {
            lNormMu <- log(x) - 0.3 ## 0.3 is HARD-CODED!
            lNormSigma <- sqrt((log(x)-lNormMu)*2)
            outputTransMatList(QmatScaled = Qmatrix, numGammaCat = 3, gammaShape = discGammaPar, numReplicates = numSamplesForTransMat, distRanFun = rlnorm, meanlog = lNormMu, sdlog = lNormSigma)
        })
        cat("Done! \n \n")
    } else{
        allIntMatList <- betweenClusTransMatList
    }
    if (missing(withinClusTransMatList)) {
        cat("Estimating transition probability matrices for branches in the within-cluster phylogenies... ")
        allExtMatList <- lapply(meanWithinBranchVec, FUN = function(x) {
            rateValue <- 1/x
            outputTransMatList(QmatScaled = Qmatrix, numGammaCat = 3, gammaShape = discGammaPar, numReplicates = numSamplesForTransMat, distRanFun = rexp, rate = rateValue)
        })
        cat("Done! \n \n")
    } else{
        allExtMatList <- betweenClusTransMatList
    }
    if (missing(poisRateNumClus)) {
        warning("No value for poisRateNumClus has been specified. Using the number of clusters in startingValues$clusInd... \n")
        poisRateNumClus <- max(startingValues$clusInd)
    } else{}

    convertedData <- .getConvertedAlignment(alignmentMat = alignment, equivVector = names(limProbs), sitePatterns = 1:ncol(alignment))
    convertedData <- lapply(convertedData, FUN = function(x) {
    colnames(x) <- rownames(alignment)
        x
    })
    seqNames <- rownames(alignment)

    logAllExtMatList <- lapply(allExtMatList, FUN = function(x) {
      lapply(x, FUN = function(y) {
        log(y)
      })
    })

    logAllIntMatList <- lapply(allIntMatList, FUN = function(x) {
      lapply(x, FUN = function(y) {
        log(y)
      })
    })

    argsForDMcore <- list(nIter = numIters, startingValues = startingValues, logLimProbs = log(limProbs), numMovesNNIint = numMovesNNIbetween, numMovesNNIext = numMovesNNIwithin, numLikThreads = numLikThreads, poisRateNumClus = poisRateNumClus, clusPhyloUpdateProp = clusPhyloUpdateProp, numSplitMergeMoves = numSplitMergeMoves, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = shiftForAlpha, logExtTransMatAll = logAllExtMatList, logIntTransMatAll = logAllIntMatList, DNAdataBin = convertedData, DNAdata = alignment)

    chainResult <- do.call(".DMphyClusCore", args = argsForDMcore)

    MAPclusInd <- chainResult[[which.max(sapply(chainResult, function(x) x$logPostProb))]]$clusInd
    list(chain = chainResult, MAPestimate = MAPclusInd)
}

#' Log-prior probability for a cluster assignment indices vector.
#'
#' This function outputs the log-prior probability for a given cluster assignment index vector, concentration parameter, and number of clusters. It is mostly used for verification purposes.
#'
#' @param clusInd vector of cluster assingment indices (should not contain any gaps)
#' @param alpha positive number, concentration parameter in the Dirichlet distribution
#' @param k positive integer giving the number of clusters; if not provided by the user, assumed to be equal to max(clusInd)
#' @return A negative numeric value.
#'
#' @examples
#' \dontrun{
#' clusIndLogPrior(clusInd = c(1,2,2,3,4,4,4), alpha = 1)
#' }
#' @export
clusIndLogPrior <- function(clusInd, alpha, k) {

    countsVec <- as.vector(table(clusInd))
    if (missing(k)) {
        k <- length(countsVec)
    } else{}
    if (length(alpha) == 1) {
        alphaVec <- rep(alpha, k)
    } else{
        alphaVec <- alpha
    }
    if (length(countsVec) < k) {
        countsVec <- c(countsVec, rep(0, k - length(countsVec)))
    } else{}
    lgamma(length(clusInd)+1) - sum(lgamma(countsVec+1)) + .logBigB(alphaVec + countsVec) - .logBigB(alphaVec)
}

#' Phylogenetic log-likelihood
#'
#' This function gives the log-likelihood for a given alignment, conditional on a phylogeny split in between- and within- cluster components.
#'
#' @param clusterPhylos list of phylo objects or NULL (indicating a singleton) with names matching tip labels in betweenPhylo
#' @param betweenPhylo phylo object with tip labels matching list names in clusterPhylos
#' @param betweenTransMatList list of matrices giving transition probabilities in the supporting phylogeny by rate variation category
#' @param withinTransMatList list of matrices giving transition probabilities in the within-cluster phylogenies by rate variation category
#' @param clusInd vector of cluster assingment indices (should not contain any gaps)
#' @param alignment a character matrix with rows corresponding to sequences and columns to sites
#' @param limProbs a vector giving limiting probabilities for all states in alignment
#' @param numLikThreads positive integer giving the number of threads used
#' @param basicChecks boolean indicating whether basic checks should be performed to ensure the call does not result in a segmentation fault
#'
#' @details This function relies on compiled C++ code and can cause segmentation faults if the supplied phylo objects are misformed. Basic checks are performed by default. The function will be faster if they are turned off.
#' @return A negative value.
#'
#' @examples
#' \dontrun{
#' INPUT_EXAMPLE()
#' }
#'
#' @export
logLikFromSplitPhylo <- function(clusterPhylos, betweenPhylo, betweenTransMatList, withinTransMatList, clusInd, alignment, limProbs, numLikThreads = 1, basicChecks = TRUE) {
    if (basicChecks) {
        environment(.checkArgumentsLogLikFromSplitPhylo) <- environment()
        .checkArgumentsLogLikFromSplitPhylo()
    } else{}
    logBetweenTransMatList <- lapply(betweenTransMatList, FUN = function(x) {
      lapply(x, FUN = function(y) {
        log(y)
      })
    })

    logWithinTransMatList <- lapply(withinTransMatList, FUN = function(x) {
      lapply(x, FUN = function(y) {
        log(y)
      })
    })

    dataBin <- .getConvertedAlignment(alignmentMat = alignment, equivVector = names(limProbs), numOpenMP = numLikThreads)
    alignmentBin <- lapply(dataBin, FUN = function(x) {
        colnames(x) <- rownames(alignment)
        x
    })
    alignmentMultiBinByClus <- lapply(names(clusterPhylos), FUN = function(x) {
        .outputDNAdataMultiBin(clusterPhylo = clusterPhylos[[x]], logLimProbs = log(limProbs), clusName = x, numLikThreads = numLikThreads, logExtMatList = logWithinTransMatList, DNAdataBin = alignmentBin, clusInd = clusInd)
    })
    names(alignmentMultiBinByClus) <- names(clusterPhylos)

    alignmentMultiBin <- redimMultiBinByClus(alignmentMultiBinByClus) ## Danger: the functions don't produce an error if we ask for alignmentMultiBin when it is not defined, because they assume that we mean alignmentMultiBinByClus.

    if (!identical(colnames(alignmentMultiBin[[1]][[1]]), betweenPhylo$tip.label)) {
        stop("The ordering of the columns in alignmentMultiBin should match that of the tip labels in betweenPhylo! Line 499.\n")
    } else{}

   .logLikCpp(edgeMat = betweenPhylo$edge, logLimProbsVec = log(limProbs), logTransMatList = logBetweenTransMatList,  numOpenMP = numLikThreads, alignmentBin = alignmentMultiBin, internalFlag = TRUE, returnRootMat = FALSE)
}

## The next function outputs a list of transition rate matrices
#' Obtaining transition probability matrices
#'
#' This function produces a Monte Carlo estimate of transition probability matrices under the assumption of discrete gamma among-sites rate variation.
#'
#' @param QmatScaled a matrix indicating transition rates among all potential states; should be scaled in such a way that each row sums to 0 and -diag(QmatScaled)*(chain limiting probabilities) = 1.
#' @param numGammaCat positive integer giving the number of rate variation categories
#' @param gammaShape positive number giving the value of the discrete gamma shape parameter
#' @param numReplicates positive integer giving the number of random values used in obtaining the estimate
#' @param distRanFun function for simulating random values, must have an argument called 'n', standing for the number of values to simulate, e.g. rnorm,rexp
#' @param ... additional parameters for distRanFun
#'
#' @details This function produces estimates of the transition probability matrices along branches in the between- or within- cluster phylogenies. It is called by DMphyClus when betweenClusTransMatList and withinClusTransMatList are not provided. The estimation can take a long time, so we recommend doing it only once with this function, saving the output, and re-using it whenever necessary.
#'
#' @return A list with length equal to numGammaCat. Each element is a matrix with dimensions equal to those of QmatScaled.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export
outputTransMatList <- function(QmatScaled, numGammaCat, gammaShape, numReplicates = 10^5, distRanFun, ...) {

    gammaRates <- phangorn::discrete.gamma(k = numGammaCat, alpha = gammaShape) ## The discrete gamma rates for discrete gamma rate variation among sites.
    lapply(gammaRates, FUN = function(gammaRate) {
        ranValues <- do.call(distRanFun, args = list(n = numReplicates, ...))
        probMatrices <- lapply(ranValues, FUN = function(x) expm::expm(QmatScaled*x*gammaRate))
        Reduce(f = "+", x = probMatrices)/length(probMatrices)
    })
}

#' Phylogenetic log-likelihood
#'
#' This function gives the log-likelihood for a given alignment, conditional on a phylogeny and a cluster assignment index vector.
#'
#' @param phylogeny phylo object
#' @param clusInd vector of cluster assingment indices (should not contain any gaps)
#' @param betweenTransMatList list of matrices giving transition probabilities in the supporting phylogeny by rate variation category
#' @param withinTransMatList list of matrices giving transition probabilities in the within-cluster phylogenies by rate variation category
#' @param alignment a character matrix with rows corresponding to sequences and columns to sites
#' @param limProbs a vector giving limiting probabilities for all states in alignment
#' @param numLikThreads positive integer giving the number of threads used
#' @param basicChecks boolean indicating whether basic checks should be performed to ensure the call does not result in a segmentation fault
#'
#' @details This function is an alternative to logLikFromSplitPhylo. Just like that function, it relies on compiled C++ code and can cause segmentation faults if the supplied phylo objects are misformed. Basic checks are performed by default. The function will be faster if they are turned off.
#' @return A negative value.
#'
#' @examples
#' \dontrun{
#' INPUT_EXAMPLE()
#' }
#' @export
logLikFromClusInd <- function(phylogeny, betweenTransMatList, withinTransMatList, clusInd, limProbs, numLikThreads, alignment, basicChecks = TRUE) {
    if (basicChecks) {
        environment(.checkArgumentsLogLikFromClusInd) <- environment()
        .checkArgumentsLogLikFromClusInd()
    } else{}

    logBetweenTransMatList <- lapply(betweenTransMatList, FUN = function(x) {
        log(x)
    })

    logWithinTransMatList <- lapply(withinTransMatList, FUN = function(x) {        log(x)
    })

    dataBin <- .getConvertedAlignment(alignmentMat = alignment, equivVector = names(limProbs), numOpenMP = numLikThreads)
    alignmentBin <- lapply(dataBin, FUN = function(x) {
        colnames(x) <- rownames(alignment)
        x
    })
    betweenPhylo <- phylogeny
    clusLabels <- as.character(1:length(unique(clusInd))) ## Labels can be numeric...
    outputClusPhySetIntPhy <- function(clusLabel, seqNames) {
        if (length(seqNames) == 1) {
            betweenPhylo$tip.label[match(seqNames, betweenPhylo$tip.label)] <<- clusLabel
            return(NA) ## The cluster is a singleton. The phylogeny for it is degenerate. The cluster root matches the tip corresponding to the sequence.
        } else{}
        clusMRCAnode <- ape::getMRCA(betweenPhylo, tip = seqNames)
        betweenPhylo$node.label[clusMRCAnode - ape::Ntip(betweenPhylo)] <<- clusLabel ## At this step, node supporting clusters are identified. All branches in clusters are given a different branch length prior. The label associates each supporting node with the subphylogenies in clusterPhylos.
        clusterPhylo <- ape::extract.clade(betweenPhylo, node = clusMRCAnode)
        ##betweenPhylo <<- ape::reorder(ape::drop.tip(betweenPhylo, seqNames, trim.internal = FALSE)) ## This leaves a tip with the cluster label. The cluster label is a number associating tips in the internal phylogeny with the phylogenies in clusterPhylos.
        newClusInd <- seq_along(betweenPhylo$tip.label)
        names(newClusInd) <- betweenPhylo$tip.label
        newClusInd[seqNames] <- 10^6
        newClusInd <- .relabel(newClusInd)
        multiForRemoval <- .introduceMultiPhyloWithDist(phylogeny = betweenPhylo, clusInd = newClusInd)
        betweenPhylo <<- ape::reorder.phylo(ape::drop.tip(multiForRemoval, seqNames, trim.internal = FALSE))
        uniquesInClus <- grepl(clusterPhylo$tip.label, pattern = "unique")
        if (any(uniquesInClus)) {
            clusterPhylo$tip.label[uniquesInClus] <- substr(clusterPhylo$tip.label[uniquesInClus], start = 1, stop = nchar(clusterPhylo$tip.label[uniquesInClus]) - 6)
        } else{}
        clusterPhylo
    }
    clusterPhylos <- mapply(clusLabels, split(names(clusInd), f = clusInd), FUN = outputClusPhySetIntPhy, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    clusterPhylos <- clusterPhylos[betweenPhylo$tip.label] ## The order of the elements in clusterPhylos must match that of the tip labels in betweenPhylo.
    oriNames <- names(clusInd)[!grepl(names(clusInd), pattern = "C")]
    names(clusInd)[!grepl(names(clusInd), pattern = "C")] <- substr(oriNames, start = 1, stop = nchar(oriNames) - 6) ## This restores the original names for clusInd. We changed them because they were responsible for confusing the function that creates betweenPhylo.

    alignmentMultiBinByClus <- lapply(names(clusterPhylos), FUN = function(x) {

        .outputDNAdataMultiBin(clusterPhylo = clusterPhylos[[x]], clusName = x, clusInd = clusInd, logExtMatList = logWithinTransMatList, numLikThreads = numLikThreads, logLimProbs = log(limProbs), DNAdataBin = alignmentBin)
    })
    names(alignmentMultiBinByClus) <- names(clusterPhylos)

    alignmentMultiBin <- redimMultiBinByClus(alignmentMultiBinByClus) ## Danger: the functions don't produce an error if we ask for alignmentMultiBin when it is not defined, because they assume that we mean alignmentMultiBinByClus.

    if (!identical(colnames(alignmentMultiBin[[1]][[1]]), betweenPhylo$tip.label)) {
        stop("The ordering of the columns in alignmentMultiBin should match that of the tip labels in betweenPhylo! Line 499.\n")
    } else{}

   .logLikCpp(edgeMat = betweenPhylo$edge, logLimProbsVec = log(limProbs), logTransMatList = logBetweenTransMatList,  numOpenMP = numLikThreads, alignmentBin = alignmentMultiBin, internalFlag = TRUE, returnRootMat = FALSE) ## Make sure this matches the result from the call to logLikCpp when the full phylogeny is used!
}
#' @useDynLib DMphyClus
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Nucleotide substitution rate matrix
#'
#' This matrix, from Posada et al. 2001, gives nucleotide substitution rates
#' in a standard HIV-1 alignment.
#'
#'
#' @format A matrix with 4 rows and 4 columns.
#' @source \url{http://mbe.oxfordjournals.org/content/18/6/897.short}
"QmatrixPosada2001"
#> [1] "QmatrixPosada2001"

#' Between-cluster transition probability matrices
#'
#' These matrices are meant to be used in the vignette example. They were
#' obtained with \code{outputTransMatList}. They represent
#' transition probabilities along branches in the between-cluster (supporting)
#' phylogeny. Among-sites rate variation follows a discrete gamma distribution
#' with 3 categories, and with parameter 0.7589 (Posada et al. 2001). We assumed
#' that branch lengths follow the log-normal distribution with means between 0.015
#' and 0.045 (10 equally spaced mean assumptions).
#'
#' @format A list of length 10. Each list element is itself a list, but with
#' length 3, corresponding to the number of rate variation categories. Each
#' element of the nested list is a 4x4 matrix giving transition probabilities
#' between the different nucleotide states.
#' @source \url{http://mbe.oxfordjournals.org/content/18/6/897.short}
"betweenClusTransProbs"
#> [1] "betweenClusTransProbs"

#' Within-cluster transition probability matrices
#'
#' These matrices are meant to be used in the vignette example.
#' They were obtained with \code{outputTransMatList}.They represent
#' transition probabilities along branches in the within-cluster phylogenies.
#' Among-sites rate variation follows a discrete gamma distribution
#' with 3 categories, and with parameter 0.7589 (Posada et al. 2001).
#' We assumed that branch lengths follow the exponential distribution
#' with means between 0.0005 and 0.0015 (10 equally spaced mean assumptions).
#'
#' @format A list of length 10. Each list element is itself a list, but with
#' length 3, corresponding to the number of rate variation categories. Each
#' element of the nested list is a 4x4 matrix giving transition probabilities
#' between the different nucleotide states.
#' @source \url{http://mbe.oxfordjournals.org/content/18/6/897.short}
"withinClusTransProbs"
#> [1] "withinClusTransProbs"

#' An alignment of six HIV-1 DNA sequences in DNAbin format
#'
#' This object comprises six HIV-1 sequences of length 1617, downloaded
#' from the Los Alamos HIV sequence database (accession numbers AB485638,
#' AB485639, AB485640, AB220944, AB220945, AB220946). Three sequences
#' are subtype B, while the three remaining sequences are CRF-AE. They
#' contain no missing data or ambiguities. We use them in the vignette.
#'
#' @format A DNAbin object.
#' @source \url{https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html}
"seqsToCluster"
#> [1] "seqsToCluster"
