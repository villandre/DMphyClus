.performStepPhylo <- function(currentValue, limProbs, shapePriorAlpha, scalePriorAlpha, numGammaCat, intTransMatAll, extTransMatAll, enableJainNeal = FALSE, currentIter, DNAdata, numNeighbours, numMovesNNIint = 1, numMovesNNIext = 1, numLikThreads, DNAdataBin, singletonMatrices, poisRateNumClus, clusPhyloUpdateProp = 1, alphaMin, numSplitMergeMoves) {
    if (class(currentValue$paraValues$internalPhylo) == "phylo") {
      currentValue <- .updateTransMatrix(currentValue = currentValue, allTransMatList = intTransMatAll, samplingRadius = 2, internalBool = TRUE, limProbs = limProbs, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin)
    }
    currentValue <- .updateTransMatrix(currentValue = currentValue, allTransMatList = extTransMatAll, currentIntTransMatList = intTransMatAll[[currentValue$paraValues$intMatListIndex]], samplingRadius = 2, internalBool = FALSE, limProbs = limProbs, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin)

    if (length(currentValue$paraValues$clusterPhylos) > 2) {
        internalPhyloList <- .updateSupportingPhylo(internalPhylo = currentValue$paraValues$internalPhylo, limProbs = limProbs, transMatList = intTransMatAll[[currentValue$paraValues$intMatListIndex]], numMovesNNI = numMovesNNIint, numLikThreads = numLikThreads, DNAdataMultiBinByClus = currentValue$DNAdataMultiBinByClus, currentLogLik = currentValue$logLik) ## This step updates the supporting phylogeny.
        currentValue$logPostProb <- currentValue$logPostProb - currentValue$logLik + internalPhyloList$logLik ## Updating the internal phylogeny doesn't change the prior value.
        currentValue$logLik <- internalPhyloList$logLik
        currentValue$paraValues$internalPhylo <- internalPhyloList$internalPhylo
    } else{}

    currentValue <- .updateClusterPhylos(currentValue = currentValue, limProbs = limProbs, intTransMatList = intTransMatAll[[currentValue$paraValues$intMatListIndex]], clusTransMatList = extTransMatAll[[currentValue$paraValues$extMatListIndex]], numMovesNNI = numMovesNNIext, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin, clusPhyloUpdateProp = clusPhyloUpdateProp) ## This step updates the cluster-specific phylogenies. It is the slowest step.
    lapply(1:numSplitMergeMoves, FUN = function(x) { ##
        currentValue <<- .splitJoinClusterMove(currentValue = currentValue, DNAdataBin = DNAdataBin, limProbs = limProbs, clusTransMatList = extTransMatAll[[currentValue$paraValues$extMatListIndex]], intTransMatList = intTransMatAll[[currentValue$paraValues$intMatListIndex]], numLikThreads = numLikThreads, singletonMatrices = singletonMatrices, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapePriorAlpha, scaleForAlpha = scalePriorAlpha, alphaMin = alphaMin)
    })

    currentValue <- .updateAlpha(currentValue, shapePriorAlpha, scalePriorAlpha, alphaMin = alphaMin)

    currentValue
}

.updateAlpha <- function(currentValue, shapeForAlpha, scaleForAlpha, alphaMin = 6) {
    increment <- (runif(1) - 0.5)
    newAlpha <- currentValue$paraValues$alpha + increment
    newAlpha <- newAlpha + abs(newAlpha - alphaMin)*(newAlpha < alphaMin)
    newLogPrior <- clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = newAlpha) + dgamma(newAlpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)
    MHratio <- exp(newLogPrior - (clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)))
    if (runif(1) < MHratio) {
        currentValue$paraValues$alpha <- newAlpha
        currentValue$logPostProb <- currentValue$logLik + newLogPrior
    } else{}
    currentValue
}

.initCurrentValue <- function(startingValues, DNAdataBin, intTransMatAll, extTransMatAll, limProbs, numLikThreads, shapeForAlpha, scaleForAlpha, alphaMin) {
    currentValue <- list()
    startingValues$phylogeny$node.label <- NULL

    if (!all(is.na(suppressWarnings(as.numeric(startingValues$phylogeny$tip.label))))) {
      warning("Not a good idea to have numbers as sequence names... The algorithm will still run (by adding 'unique' to their names to prevent confusion), but change that! \n")
      numberPosPhylo <- which(!is.na(suppressWarnings(as.numeric(startingValues$phylogeny$tip.label))))
      numberPosClusInd<- which(!is.na(suppressWarnings(as.numeric(names(startingValues$clusInd)))))
      startingValues$phylogeny$tip.label[numberPosPhylo] <- paste(startingValues$phylogeny$tip.label[numberPosPhylo], "unique", sep = "")
      startingValues$clusInd[numberPosClusInd] <- paste(startingValues$clusInd[numberPosClusInd], "unique", sep = "")
    } else{}
    currentValue$paraValues <- startingValues ## Initializing starting values. Note that the cluster indices are initialized in the wrapper function onePolyaIter.
    currentValue$paraValues$intMatListIndex <- ceiling(length(intTransMatAll)/2)
    currentValue$paraValues$extMatListIndex <- ceiling(length(extTransMatAll)/2)
    currentValue$paraValues$clusInd <- startingValues$clusInd[startingValues$phylogeny$tip.label]

    # names(currentValue$paraValues$clusInd)[!grepl(names(currentValue$paraValues$clusInd), pattern = "C")] <- paste(names(currentValue$paraValues$clusInd)[!grepl(names(currentValue$paraValues$clusInd), pattern = "C")], "unique", sep = "")
    if (max(currentValue$paraValues$clusInd) == 1) {
      internalPhylo <- NULL
      currentValue$paraValues$clusterPhylos[[1]] <- currentValue$paraValues$phylogeny
      names(currentValue$paraValues$clusterPhylos) <- "1"
      DNAdataMultiBinByClus <- list(.outputDNAdataMultiBin(clusterPhylo = currentValue$paraValues$clusterPhylos[[1]], clusName = "1", clusInd = currentValue$paraValues$clusInd, DNAdataBin = lapply(DNAdataBin, FUN = function(x) x[,currentValue$paraValues$clusterPhylos[[1]]$tip.label]), extMatList = extTransMatAll[[currentValue$paraValues$extMatListIndex]], numLikThreads = numLikThreads, limProbs = limProbs))
      names(DNAdataMultiBinByClus) <- names(currentValue$paraValues$clusterPhylos)
      currentValue$DNAdataMultiBinByClus <- DNAdataMultiBinByClus ## Will need to be adjusted with redimMultiBinByClus before being used in a call to logLikCpp
      currentValue$logLik <- .logLikCpp(edgeMat = currentValue$paraValues$clusterPhylos[[1]]$edge, limProbsVec = limProbs, transMatList = extTransMatAll[[currentValue$paraValues$extMatListIndex]],  numOpenMP = numLikThreads, alignmentBin = lapply(DNAdataBin, FUN = function(x) {x[,currentValue$paraValues$clusterPhylos[[1]]$tip.label]}), internalFlag = FALSE, returnRootMat = FALSE)
    } else {
      internalPhylo <- currentValue$paraValues$phylogeny
      clusLabels <- as.character(1:length(unique(currentValue$paraValues$clusInd)))
      outputClusPhySetIntPhy <- function(clusLabel, seqNames) {
        if (length(seqNames) == 1) {
          internalPhylo$tip.label[match(seqNames, internalPhylo$tip.label)] <<- clusLabel
          return(NA) ## The cluster is a singleton. The phylogeny for it is degenerate. The cluster root matches the tip corresponding to the sequence.
        } else{}
        clusMRCAnode <- ape::getMRCA(internalPhylo, tip = seqNames)
        internalPhylo$node.label[clusMRCAnode - ape::Ntip(internalPhylo)] <<- clusLabel ## At this step, node supporting clusters are identified. All branches in clusters are given a different branch length prior. The label associates each supporting node with the subphylogenies in clusterPhylos.
        clusterPhylo <-ape::extract.clade(internalPhylo, node = clusMRCAnode)
        ##internalPhylo <<- ape::reorder.phylo(ape::drop.tip(internalPhylo, seqNames, trim.internal = FALSE)) ## This leaves a tip with the cluster label. The cluster label is a number associating tips in the internal phylogeny with the phylogenies in clusterPhylos.
        newClusInd <- seq_along(internalPhylo$tip.label)
        names(newClusInd) <- internalPhylo$tip.label
        newClusInd[seqNames] <- 10^6
        newClusInd <- .relabel(newClusInd)
        multiForRemoval <- .introduceMultiPhyloWithDist(phylogeny = internalPhylo, clusInd = newClusInd)
        internalPhylo <<- ape::reorder.phylo(ape::drop.tip(multiForRemoval, seqNames, trim.internal = FALSE))
        uniquesInClus <- grepl(clusterPhylo$tip.label, pattern = "unique")
        if (any(uniquesInClus)) {
          clusterPhylo$tip.label[uniquesInClus] <- substr(clusterPhylo$tip.label[uniquesInClus], start = 1, stop = nchar(clusterPhylo$tip.label[uniquesInClus]) - 6)
        } else{}
        clusterPhylo
      }
      currentValue$paraValues$clusterPhylos <- mapply(clusLabels, split(names(currentValue$paraValues$clusInd), f = currentValue$paraValues$clusInd), FUN = outputClusPhySetIntPhy, SIMPLIFY = FALSE, USE.NAMES = TRUE)
      currentValue$paraValues$clusterPhylos <- currentValue$paraValues$clusterPhylos[internalPhylo$tip.label] ## The order of the elements in clusterPhylos must match that of the tip labels in internalPhylo.
      currentValue$paraValues$internalPhylo <- internalPhylo

      # oriNames <- names(currentValue$paraValues$clusInd)[!grepl(names(currentValue$paraValues$clusInd), pattern = "C")]
      # names(currentValue$paraValues$clusInd)[!grepl(names(currentValue$paraValues$clusInd), pattern = "C")] <- substr(oriNames, start = 1, stop = nchar(oriNames) - 6) ## This restores the original names for clusInd. We changed them because they were responsible for confusing the function that creates internalPhylo.

      DNAdataMultiBinByClus <- lapply(names(currentValue$paraValues$clusterPhylos), FUN = function(x) {
        .outputDNAdataMultiBin(clusterPhylo = currentValue$paraValues$clusterPhylos[[x]], clusName = x, clusInd = currentValue$paraValues$clusInd, DNAdataBin = DNAdataBin, extMatList = extTransMatAll[[currentValue$paraValues$extMatListIndex]], numLikThreads = numLikThreads, limProbs = limProbs)
      })
      names(DNAdataMultiBinByClus) <- names(currentValue$paraValues$clusterPhylos)

      currentValue$DNAdataMultiBinByClus <- DNAdataMultiBinByClus ## Will need to be adjusted with redimMultiBinByClus before being used in a call to logLikCpp
      DNAdataMultiBin <- redimMultiBinByClus(currentValue$DNAdataMultiBinByClus) ## Danger: the functions don't produce an error if we ask for DNAdataMultiBin when it is not defined, because they assume that we mean DNAdataMultiBinByClus.
      if (!identical(colnames(DNAdataMultiBin[[1]][[1]]), currentValue$paraValues$internalPhylo$tip.label)) {
        stop("The ordering of the columns in DNAdataMultiBin should match that of the tip labels in internalPhylo! Line 499.\n")
      } else{}
      currentValue$logLik <- .logLikCpp(edgeMat = currentValue$paraValues$internalPhylo$edge, limProbsVec = limProbs, transMatList = intTransMatAll[[currentValue$paraValues$intMatListIndex]],  numOpenMP = numLikThreads, alignmentBin = DNAdataMultiBin, internalFlag = TRUE, returnRootMat = FALSE) ## Make sure this matches the result from the call to logLikCpp when the full phylogeny is used!
    }
    currentValue$clusterCounts <- as.vector(table(currentValue$paraValues$clusInd)) ## The as.vector ensures that clusterCounts behaves always as a vector, but it removes the names.
    names(currentValue$clusterCounts) <- 1:max(currentValue$paraValues$clusInd) ## This once again rests on the assumption that there is no gap in the cluster labels. This is an assumption we made before. Names are needed sometimes, although it is better to index by number, rather than by name (must require looking into an index).

    currentValue$paraValues <- currentValue$paraValues[-match("phylogeny", names(currentValue$paraValues))] ## We remove the phylogeny element, that's now been broken into internal and external phylogenies to reduce the amount of bookkeeping required.
    currentValue$logPostProb <- currentValue$logLik + clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE) ## The starting value for log-posterior probability for the starting partition.
    currentValue
}

.performSplit <- function(currentValue, numSplits, clusLabel, DNAdataBin, limProbs, clusTransMatList, intTransMatList, numLikThreads, singletonMatrices, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    originalPhylo <- currentValue$paraValues$clusterPhylos[[clusLabel]]
    newClusLabel <- as.character(max(currentValue$paraValues$clusInd)+1)
    graftedTree <- ape::rtree(2)
    graftedTree$edge.length <- NULL
    graftedTree$tip.label <- c(clusLabel, newClusLabel)
    if (length(currentValue$paraValues$clusterPhylos) > 1) {
        newInternalTree <- ape::bind.tree(x = ape::reorder.phylo(currentValue$paraValues$internalPhylo), y = graftedTree, where = match(clusLabel, currentValue$paraValues$internalPhylo$tip.label)) ## Breaking up a cluster results in the internal phylogeny changing too, since it must now integrate the part of the cluster phylogeny that used to link the two new clusters. The ordering of the labels in the size 2 phylogeny supporting the clusters doesn't matter: it is symmetrical.
    } else {
        newInternalTree <- graftedTree
    }
    DNAdataMultiBinByClusCopy <- currentValue$DNAdataMultiBinByClus
    rootChildren <- phangorn::Children(originalPhylo, ape::Ntip(originalPhylo)+1)
    splitDNAbyClus <- lapply(rootChildren, FUN = function(childNode) {
        if (childNode <= ape::Ntip(originalPhylo)) {
            return(list(tree = NA, clusMatrix = singletonMatrices[[originalPhylo$tip.label[childNode]]]))
        } else{
            newSubTree <- ape::extract.clade(phy = originalPhylo, node = childNode)
            alignmentBin <- .dataBinSubset(DNAdataBin = DNAdataBin, keepSeqNamesOrNums = newSubTree$tip.label)

            if (!identical(colnames(alignmentBin[[1]]), newSubTree$tip.label)) {
                stop("Look at names... \n")
            } else{}
            return(list(tree = newSubTree, clusMatrix = .logLikCpp(edgeMat = newSubTree$edge, limProbsVec = limProbs, transMatList = clusTransMatList, numOpenMP = numLikThreads, alignmentBin = alignmentBin, internalFlag = FALSE, returnRootMat = TRUE)))
        }
    })
    DNAdataMultiBinByClusCopy[graftedTree$tip.label] <- lapply(splitDNAbyClus, FUN = function(x) x$clusMatrix)
    newClusterPhylo <- lapply(splitDNAbyClus, FUN = function(x) x$tree)
    newClusPhylosLabel <- graftedTree$tip.label

    DNAdataMultiBinByClusCopy <- DNAdataMultiBinByClusCopy[newInternalTree$tip.label]
    newDNAdataMultiBin <- redimMultiBinByClus(DNAdataMultiBinByClusCopy) ## This works only if the number of elements in the list on the LHS matches that in the list on the RHS.
    if (!identical(colnames(newDNAdataMultiBin[[1]][[1]]), newInternalTree$tip.label)) {
        stop("Look at column names! \n")
    } else{}
    newLogLik <- .logLikCpp(edgeMat = newInternalTree$edge, limProbsVec = limProbs, transMatList = intTransMatList, numOpenMP = numLikThreads, alignmentBin = newDNAdataMultiBin, internalFlag = TRUE, returnRootMat = FALSE)
    newClusInd <- currentValue$paraValues$clusInd
    if (!(class(splitDNAbyClus[[2]]$tree) == "phylo"))  {
        labelsToAmend <- originalPhylo$tip.label[rootChildren[[2]]]
    } else{
        labelsToAmend <- splitDNAbyClus[[2]]$tree$tip.label
    }
    newClusInd[labelsToAmend] <- max(currentValue$paraValues$clusInd) + 1 ## Sequences in the cluster supported by remainingTree are assigned a label equal to max(clusInd) + 1, while those supported by newClusTree keep the label of the cluster that was split.
    newCounts <- table(newClusInd)
    newLogPostProb <- newLogLik + clusIndLogPrior(clusInd = newClusInd, alpha = currentValue$paraValues$alpha) + dpois(length(newCounts), lambda = poisRateNumClus, log = TRUE) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)## Added the Poisson log-prob. to reflect a Poisson prior on the total number of clusters.

    tipAncestorsNew <- phangorn::Ancestors(newInternalTree, node = 1:ape::Ntip(newInternalTree), type = "parent")
    names(tipAncestorsNew) <- newInternalTree$tip.label
    ancestorsGroupsNew <- split(tipAncestorsNew, f = tipAncestorsNew)
    ancestorsGroupsPairsNew <- ancestorsGroupsNew[sapply(ancestorsGroupsNew, FUN = function(x) {length(x) > 1})]
    numPairsNew <- length(ancestorsGroupsPairsNew)
    transKernRatio <- numSplits/numPairsNew
    list(logLik = newLogLik, clusInd = newClusInd, counts = newCounts, logPostProb = newLogPostProb, transKernRatio = transKernRatio, clusterPhylos = newClusterPhylo, internalPhylo = newInternalTree, DNAdataMultiBinByClus = DNAdataMultiBinByClusCopy, newClusPhylosLabel = newClusPhylosLabel)
}

.performMerge <- function(currentValue, ancestorsGroupsPairs, pairToSelect, DNAdataBin, limProbs, clusTransMatList, intTransMatList, numLikThreads, singletonMatrices, tipAncestors, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    numPairs <- length(ancestorsGroupsPairs)
    clustersInPair <- names(ancestorsGroupsPairs[[pairToSelect]])
    hostTree <- ape::rtree(2, tip.label = c(0,0))
    hostTree$edge.length <- 0
    if (class(currentValue$paraValues$clusterPhylos[[clustersInPair[[1]]]]) == "phylo") {
        hostTree <- ape::bind.tree(x = hostTree, y = ape::reorder.phylo(currentValue$paraValues$clusterPhylos[[clustersInPair[[1]]]]), where = 1)
        hostChildren <- phangorn::Children(hostTree, ape::Ntip(hostTree)+1)
        whereToBind <- hostChildren[hostChildren < ape::Ntip(hostTree)]
    } else{
        hostTree$tip.label[1] <- names(currentValue$paraValues$clusInd)[match(as.numeric(clustersInPair[[1]]), currentValue$paraValues$clusInd)]
        whereToBind <- 2
    }

    if (class(currentValue$paraValues$clusterPhylos[[clustersInPair[[2]]]]) == "phylo") {
        newClusterPhylo <- ape::bind.tree(x = hostTree, y = ape::reorder.phylo(currentValue$paraValues$clusterPhylos[[clustersInPair[[2]]]]), where = whereToBind)
    } else{
        hostTree$tip.label[whereToBind] <- names(currentValue$paraValues$clusInd)[match(as.numeric(clustersInPair[[2]]), currentValue$paraValues$clusInd)]
        newClusterPhylo <- hostTree
    }
    newTipLabel <- as.character(max(currentValue$paraValues$clusInd)+1)
    currentValue$paraValues$internalPhylo$node.label[[tipAncestors[[clustersInPair[[1]]]]-ape::Ntip(currentValue$paraValues$internalPhylo)]] <- newTipLabel ## The new cluster will be given the label of a totally new cluster.
    newInternalTree <- suppressWarnings(ape::drop.tip(currentValue$paraValues$internalPhylo, tip = clustersInPair, trim.internal = FALSE)) ## The node label should become a tip label. There will be a gap in the nodes. We suppress warnings, because the algorithm can handle the case where all tips are dropped.
    mergeDataBin <- .dataBinSubset(DNAdataBin = DNAdataBin, keepSeqNamesOrNums = newClusterPhylo$tip.label)
    DNAdataMultiBinByClusCopy <- currentValue$DNAdataMultiBinByClus[-match(clustersInPair, names(currentValue$DNAdataMultiBinByClus))]
    clusTransMatSelectIndex <- min(length(clusTransMatList), ape::Ntip(newClusterPhylo))
    if (!identical(newClusterPhylo$tip.label, colnames(mergeDataBin[[1]]))) {
        stop("Look at colnames! \n")
    } else{}
    DNAdataMultiBinByClusCopy[[newTipLabel]] <- .logLikCpp(edgeMat = newClusterPhylo$edge, limProbsVec = limProbs, transMatList = clusTransMatList, numOpenMP = numLikThreads, alignmentBin = mergeDataBin, internalFlag = FALSE, returnRootMat = TRUE) ## This is automatically added to the end of the list, which means the ordering of the list should not be affected.

    newClusPhylosLabel <- newTipLabel
    newDNAdataMultiBin <- redimMultiBinByClus(DNAdataMultiBinByClusCopy)
    if (!is.null(newInternalTree)) {
      if (!identical(newInternalTree$tip.label, colnames(newDNAdataMultiBin[[1]][[1]]))) {
        stop("Look at colnames! \n")
      } else{}
      newLogLik <- .logLikCpp(edgeMat = newInternalTree$edge, limProbsVec = limProbs, transMatList = intTransMatList, numOpenMP = numLikThreads, alignmentBin = newDNAdataMultiBin, internalFlag = TRUE, returnRootMat = FALSE)
    } else {
      if (!identical(colnames(mergeDataBin[[1]]), newClusterPhylo$tip.label)) {
        stop("Error in .performMerge: orders of columns in mergeDataBin should match that of the tip labels in newClusterPhylo!\n")
      } else{}
      newLogLik <- .logLikCpp(edgeMat = newClusterPhylo$edge, limProbsVec = limProbs, transMatList = clusTransMatList, numOpenMP = numLikThreads, alignmentBin = mergeDataBin, internalFlag = FALSE, returnRootMat = FALSE) ## Should make sure that the order in mergeDataBin matches that in newClusterPhylo...
    }
    newClusInd <- replace(currentValue$paraValues$clusInd, which(currentValue$paraValues$clusInd %in% as.numeric(clustersInPair)), as.numeric(newTipLabel)) ## There's a gap in the cluster labels.
    newCounts <- table(newClusInd)
    newLogPostProb <- newLogLik + clusIndLogPrior(clusInd = newClusInd, alpha = currentValue$paraValues$alpha) + dpois(length(newCounts), lambda = poisRateNumClus, log = TRUE) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)## Added the Poisson log-prob. to reflect a Poisson prior on the total number of clusters.
    transKernRatio <- numPairs/sum(newCounts>1)
    newClusterPhylo <- list(newClusterPhylo)
    list(logLik = newLogLik, clusInd = newClusInd, counts = newCounts, logPostProb = newLogPostProb, transKernRatio = transKernRatio, clusterPhylos = newClusterPhylo, internalPhylo = newInternalTree, DNAdataMultiBinByClus = DNAdataMultiBinByClusCopy, newClusPhylosLabel = newClusPhylosLabel)
}

.splitJoinClusterMove <- function(currentValue, DNAdataBin, limProbs, clusTransMatList, intTransMatList, numLikThreads, singletonMatrices, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    if (class(currentValue$paraValues$internalPhylo) == "phylo") {
      tipAncestors <- phangorn::Ancestors(currentValue$paraValues$internalPhylo, node = seq_along(currentValue$paraValues$internalPhylo$tip.label), type = "parent")
      names(tipAncestors) <- currentValue$paraValues$internalPhylo$tip.label
      ancestorsGroups <- split(tipAncestors, f = tipAncestors)
      ancestorsGroupsPairs <- ancestorsGroups[sapply(ancestorsGroups, FUN = function(x) {length(x) > 1})]
      numPairs <- length(ancestorsGroupsPairs)
    } else {
      numPairs <- 0
    }
    clustersToSplit <- names(currentValue$clusterCounts)[currentValue$clusterCounts > 1]

    if ((numPairs > 0) & (length(clustersToSplit) > 0)) {
        splitMove <- runif(1)<0.5
    } else if (numPairs > 0) {
        splitMove <- FALSE
    } else if (length(clustersToSplit) > 0) {
        splitMove <- TRUE
    } else {
        return(currentValue) ## No move is possible in this setting, so the algorithm simply omits this step.
    }

    if (splitMove) {
        clusLabel <- sample(clustersToSplit, size = 1)
        splitMergeResult <- .performSplit(currentValue = currentValue, numSplits = length(clustersToSplit), clusLabel = clusLabel, DNAdataBin = DNAdataBin, limProbs = limProbs, clusTransMatList = clusTransMatList, intTransMatList = intTransMatList, numLikThreads = numLikThreads, singletonMatrices = singletonMatrices, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
    } else {
        pairToSelect <- sample.int(n = numPairs, size = 1)
        clustersInPair <- names(ancestorsGroupsPairs[[pairToSelect]])
        splitMergeResult <- .performMerge(currentValue = currentValue, ancestorsGroupsPairs = ancestorsGroupsPairs, pairToSelect = pairToSelect, DNAdataBin = DNAdataBin, limProbs = limProbs, clusTransMatList = clusTransMatList, intTransMatList = intTransMatList, numLikThreads = numLikThreads, singletonMatrices = singletonMatrices, tipAncestors = tipAncestors, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
    }

    MHratio <- splitMergeResult$transKernRatio*exp(splitMergeResult$logPostProb - currentValue$logPostProb)

    testUnif <- runif(1)
    if (testUnif < MHratio) {
        currentValue$logLik <- splitMergeResult$logLik
        currentValue$logPostProb <- splitMergeResult$logPostProb
        currentValue$DNAdataMultiBinByClus <- splitMergeResult$DNAdataMultiBinByClus
        currentValue$paraValues$clusInd <- splitMergeResult$clusInd
        currentValue$clusterCounts <- splitMergeResult$counts
        currentValue$paraValues$internalPhylo <- splitMergeResult$internalPhylo
        currentValue$paraValues$clusterPhylos[splitMergeResult$newClusPhylosLabel] <- splitMergeResult$clusterPhylos

        if (!splitMove) {
            currentValue$paraValues$clusterPhylos <- currentValue$paraValues$clusterPhylos[-match(clustersInPair, names(currentValue$paraValues$clusterPhylos))]
            names(currentValue$DNAdataMultiBinByClus) <- as.character(.relabel(as.numeric(names(currentValue$DNAdataMultiBinByClus))))
            names(currentValue$paraValues$clusterPhylos) <- as.character(.relabel(as.numeric(names(currentValue$paraValues$clusterPhylos))))
            names(currentValue$clusterCounts) <- as.character(.relabel(as.numeric(names(currentValue$clusterCounts))))
            currentValue$paraValues$clusInd <- .relabel(currentValue$paraValues$clusInd)
            currentValue$paraValues$internalPhylo$tip.label <- as.character(.relabel(as.numeric(currentValue$paraValues$internalPhylo$tip.label)))
        } else {}
    } else{}
    currentValue
}

.updateTransMatrix <- function(currentValue, allTransMatList, currentIntTransMatList, samplingRadius = 2, internalBool, limProbs, numLikThreads, DNAdataBin) {
    if (internalBool) {
        basicIndex <- currentValue$paraValues$intMatListIndex
        DNAdataMultiBin <- redimMultiBinByClus(currentValue$DNAdataMultiBinByClus)
    } else{
        basicIndex <- currentValue$paraValues$extMatListIndex
        intTransMatList <- currentIntTransMatList
    }
    newIndex <- basicIndex + sample(c(-samplingRadius:-1,1:samplingRadius), size = 1)
    newIndex <- (newIndex<=0)*(length(allTransMatList) + newIndex) + (newIndex > length(allTransMatList))*(newIndex - length(allTransMatList)) + ((newIndex > 0) & (newIndex <= length(allTransMatList)))*newIndex  ## We have a circular transition kernel. If the range for allTransMatList is wide enough, the circularity should be never invoked.
    if (internalBool) {
        intTransMatList <- allTransMatList[[newIndex]]
    } else{ ## THis is for the external phylogenies: we need to rederive DNAdataMultiBinByClus

        DNAdataMultiBinByClus <- lapply(names(currentValue$paraValues$clusterPhylos), FUN = function(x) {
            .outputDNAdataMultiBin(clusterPhylo = currentValue$paraValues$clusterPhylos[[x]], clusInd = currentValue$paraValues$clusInd, clusName = x, DNAdataBin = DNAdataBin, numLikThreads = numLikThreads, extMatList = allTransMatList[[newIndex]], limProbs = limProbs)
        })
        DNAdataMultiBin <- redimMultiBinByClus(DNAdataMultiBinByClus)
    }
    if (class(currentValue$paraValues$internalPhylo) == "phylo") {
      newLogLik <- .logLikCpp(edgeMat = currentValue$paraValues$internalPhylo$edge, limProbsVec = limProbs, transMatList = intTransMatList, numOpenMP = numLikThreads, alignmentBin = DNAdataMultiBin, internalFlag = TRUE)
    } else {
      newLogLik <- .logLikCpp(edgeMat = currentValue$paraValues$clusterPhylos[[1]]$edge, limProbsVec = limProbs, transMatList = allTransMatList[[newIndex]], numOpenMP = numLikThreads, alignmentBin = lapply(DNAdataBin, FUN = function(x) {x[,currentValue$paraValues$clusterPhylos[[1]]$tip.label]}), internalFlag = FALSE) ## Not very efficient...
    }
    MHratio <- exp(newLogLik - currentValue$logLik)
    if (runif(1) < MHratio) {
        if (internalBool) {
            currentValue$paraValues$intMatListIndex <- newIndex
        } else{
            currentValue$paraValues$extMatListIndex <- newIndex
        }
        currentValue$logPostProb <- newLogLik +  (currentValue$logPostProb - currentValue$logLik)
        currentValue$logLik <- newLogLik
    } else{}
    currentValue
}

.outputDNAdataMultiBin <- function(clusterPhylo, clusInd, clusName, DNAdataBin, extMatList, numLikThreads, limProbs) {
    if (class(clusterPhylo) != "phylo") {
        tipInClus <- names(clusInd)[match(as.numeric(clusName), clusInd)]
        alignRedim <- sapply(1:length(DNAdataBin), FUN = function(locusIndex) {
            DNAdataBin[[locusIndex]][,tipInClus]*1e300 ## The factor prevents computational zeros.
        })
        return(replicate(n = length(extMatList), expr = identity(alignRedim), simplify = FALSE))
    } else{}
    subDNAdataBin <- .dataBinSubset(DNAdataBin = DNAdataBin, keepSeqNamesOrNums = clusterPhylo$tip.label) ## Careful: the column order in the matrices of subDNAdataBin must match the ordering of the tip labels in clusterPhylos[[x]]: tested, they seem to match.
    if (!identical(colnames(subDNAdataBin[[1]]), clusterPhylo$tip.label)) {
        stop("The ordering of the columns in DNAdataMultiBin should match that of the tip labels in internalPhylo! Line 486.\n")
    } else{}

    .logLikCpp(edgeMat = clusterPhylo$edge, limProbsVec = limProbs, transMatList = extMatList, numOpenMP = numLikThreads, alignmentBin = subDNAdataBin, internalFlag = FALSE, returnRootMat = TRUE) ## The cluster specific phylogenies all have the same branch length priors across all branches.
}

.DMphyClusCore <- function(nIter, startingValues, DNAdata, limProbs, shapeForAlpha, scaleForAlpha, numMovesNNIint, numMovesNNIext, numLikThreads, DNAdataBin, poisRateNumClus, clusPhyloUpdateProp, numSplitMergeMoves, alphaMin, extTransMatAll, intTransMatAll) {

    currentValue <- .initCurrentValue(startingValues = startingValues, DNAdataBin = DNAdataBin, extTransMatAll = extTransMatAll, intTransMatAll = intTransMatAll, limProbs = limProbs, numLikThreads = numLikThreads, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
    if (is.matrix(extTransMatAll[[1]])) {
        numGammaCat <- length(extTransMatAll) ## We only have one set of substitution rate matrices.
    } else{
        numGammaCat <- length(extTransMatAll[[1]])
    }
    singletonMatrices <- lapply(names(currentValue$paraValues$clusInd), FUN = function(obsName) {
        alignRedim <- sapply(seq_along(DNAdataBin), FUN = function(locusIndex) {
            DNAdataBin[[locusIndex]][,obsName]*1e300 ## This is an offset to prevent computational zeros.
        })
        rep(list(alignRedim), numGammaCat)
    })
    names(singletonMatrices) <- names(currentValue$paraValues$clusInd)
    cat("Launching chain... \n \n")

    progress <- txtProgressBar(style = 3)
    setTxtProgressBar(pb = progress, value = 0)
    onePerc <- nIter/100
    longOut <- lapply(1:nIter, FUN = function(x) {
        if ((x %% onePerc) == 0) {
            setTxtProgressBar(pb = progress, value = x/nIter)
        } else{}

        currentValue <<- .performStepPhylo(currentValue = currentValue, limProbs = limProbs, shapePriorAlpha = shapeForAlpha, scalePriorAlpha = scaleForAlpha, extTransMatAll = extTransMatAll, intTransMatAll = intTransMatAll, currentIter = x, numMovesNNIint = numMovesNNIint, numMovesNNIext = numMovesNNIext, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin, singletonMatrices = singletonMatrices, poisRateNumClus = poisRateNumClus, clusPhyloUpdateProp = clusPhyloUpdateProp, alphaMin = alphaMin, numSplitMergeMoves = numSplitMergeMoves)

        paraVec <- currentValue$paraValues
        names(paraVec) <- replace(names(paraVec), match(c("internalPhylo", "extMatListIndex","intMatListIndex"), names(paraVec)), c("supportPhylo", "withinTransMatListIndex", "betweenTransMatListIndex"))
        c(paraVec, list(logPostProb = currentValue$logPostProb), list(logLik = currentValue$logLik))
    })
    setTxtProgressBar(pb = progress, value = 1)
    close(con = progress)
    cat("\n Chain complete. \n\n\n")


    longOut
}

.relabel <- function(xVector) {
  orderElements <- order(order(xVector))
  sortedElements <- sort(xVector)
  relabelled <- cumsum(diff(c(0,sortedElements))>0)
  relabelled[orderElements] ## Should restore the original order, but with the new labels.
}

.logBigB <- function(x) {
    logNumer <- sum(lgamma(x))
    logDenum <- lgamma(sum(x))
    logNumer - logDenum
}

## Note: If alignmentMat is a nucleotide alignment with ambiguities, standard notation must be used, i.e. nucleotides are given by  "a", "t", "c", "g", and ambiguities, by "r", "y", etc.. All codes are in small letters.
## numStates is 4 for nucleotide alignments and 20 for AA alignments.
## childNodeInClusIndic is by default 1000, meaning that branches leading to tips, and only those branches, belong to clusters.
## NEXT UPDATE: The first argument should not be edgeMat: it's not straightforward. It should be the phylo object itself. It could be passed as a Rcpp::List to C++. We could then remove the NnodeCons argument. We could also add to the matrix an element named childNodeInClusIndic, which would remove the need for the user to specify it in a separate argument.
.logLikCpp <- function(edgeMat, limProbsVec, alignmentMat, transMatList, numOpenMP, equivVector, alignmentBin, childNodeInClusIndic = 1000, returnRootMat = FALSE, internalFlag = FALSE) {
    ## The placeholders exist because Rcpp will want these arguments to be defined.
    numStatesCons <- length(limProbsVec)
    placeholderAlignmentNum <- matrix(0,1,1)
    placeholderAlignmentAlpha <- matrix("a",1,1)
    placeholderBin <- list(matrix(1000,1,1))
    placeholderEquiv <- rep("a",numStatesCons)
    if (missing(alignmentBin)) {
        if (is.character(alignmentMat[1,1])) {## Alphanumeric alignment, the function will need to convert it. It is fed to alignmentAlphaMat in logLikCppToWrap.
            output <- logLikCppToWrap(edgeMat = edgeMat, alignmentMat = placeholderAlignmentNum, limProbsVec = limProbsVec, transMatList = transMatList, numOpenMP = numOpenMP, equivVector = equivVector, alignmentAlphaMat = alignmentMat, alignmentBin = placeholderBin, childNodeInClusIndic = childNodeInClusIndic, returnMatIndic = returnRootMat, internalFlag = internalFlag)
        } else {
            output <- logLikCppToWrap(edgeMat = edgeMat, alignmentMat = alignmentMat, limProbsVec = limProbsVec, transMatList = transMatList, numOpenMP = numOpenMP, equivVector = placeholderEquiv, alignmentAlphaMat = placeholderAlignmentAlpha, alignmentBin = placeholderBin, childNodeInClusIndic = childNodeInClusIndic, returnMatIndic = returnRootMat, internalFlag = internalFlag) ## Both equivVector and alignmentAlphaMat are given arbitrary values serving as placeholders, since the underlying C++ function cannot handle missing arguments.
        }
    } else {
        output <- logLikCppToWrap(edgeMat = edgeMat, alignmentMat = placeholderAlignmentNum, limProbsVec = limProbsVec, transMatList = transMatList, numOpenMP = numOpenMP, equivVector = placeholderEquiv, alignmentAlphaMat = placeholderAlignmentAlpha, alignmentBin = alignmentBin, childNodeInClusIndic = childNodeInClusIndic, returnMatIndic = returnRootMat, internalFlag = internalFlag)
    }
    output
}

.getConvertedAlignment <- function(alignmentMat, numStatesCons, equivVector, numOpenMP = 2) {

    placeholderAllIntegers <- 2
    placeholderAlignmentNum <- matrix(0,1,1)
    placeholderAlignmentAlpha <- matrix("a",1,1)
    placeholderEquiv <- 1:numStatesCons

    if (is.character(alignmentMat[1,1])) {
       output <- getConvertedAlignmentToWrap(alignmentMat = placeholderAlignmentNum, numStatesCons = numStatesCons, numOpenMP = numOpenMP, equivVector = equivVector, alignmentAlphaMat = alignmentMat)
    } else {
        output <- getConvertedAlignmentToWrap(alignmentMat = alignmentMat, numStatesCons = numStatesCons, numOpenMP = numOpenMP, equivVector = placeholderEquiv, alignmentAlphaMat = placeholderAlignmentAlpha)  ## Both equivVector and alignmentAlphaMat are given arbitrary values serving as placeholders, since the underlying C++ function cannot handle missing arguments.
    }
    output
}

.introduceMultiPhyloWithDist <-  function(phylogeny, clusInd) {
    if (is.null(phylogeny$tip.label)) {
        toSplit <- as.character(1:ape::Ntip(phylogeny))
    } else {
        toSplit <- phylogeny$tip.label
    }
    nodeListCut <- split(x = toSplit, f = clusInd)
    keepIndices <- sapply(nodeListCut, FUN = function(x) length(x) > 2)
    nodeListCut <- nodeListCut[keepIndices]
    newPhylo <- phylogeny
    internalFun <- function(seqsInClusInd) {
        multiTips <- which(newPhylo$tip.label %in% seqsInClusInd)
        mrcaNode <- ape::getMRCA(phy = newPhylo, tip = multiTips)
        distsToMRCA <- as.vector(ape::dist.nodes(newPhylo)[mrcaNode, match(seqsInClusInd, newPhylo$tip.label)])
        newPhylo$tip.label[multiTips] <- "removeMe"
        newTree <- ape::stree(length(seqsInClusInd), tip.label = seqsInClusInd)
        newTree$edge.length <- distsToMRCA
        incrementedTree <- ape::reorder.phylo(ape::bind.tree(x = newPhylo, y = newTree, where = mrcaNode))
        newPhylo <<- ape::reorder.phylo(ape::drop.tip(phy = incrementedTree, tip = which(incrementedTree$tip.label == "removeMe")))
    }
    lapply(nodeListCut, FUN = internalFun)
    ape::reorder.phylo(newPhylo)
}

.dataBinSubset <- function(DNAdataBin, keepLociNums, keepSeqNamesOrNums) {
    DNAdataBinSub <- DNAdataBin
    if (!missing(keepLociNums)) {
        DNAdataBinSub <- DNAdataBin[keepLociNums]
    } else{}
    if (!missing(keepSeqNamesOrNums)) {
        if (class(keepSeqNamesOrNums) == "character") {
            if (is.null(colnames(DNAdataBin[[1]]))) {
                stop("Matrices in DNAdataBin should have named columns. \n")
            } else {
                keepVec <- match(keepSeqNamesOrNums, colnames(DNAdataBinSub[[1]]))
                DNAdataBinSub <- lapply(DNAdataBinSub, FUN = function(x) {
                    x[,keepVec]
                })
            }
        } else{}
    }
    DNAdataBinSub
}

## DNAdataMultiBin is a list of lists of matrices. Outer list has # elements = # rate categories. The next level has # elements = #loci, the inner level is a matrix with # rows = # states and # col. = number of tips.
.updateSupportingPhylo <- function(internalPhylo, currentLogLik, limProbs, transMatList, numMovesNNI, numLikThreads, DNAdataMultiBinByClus) {

    DNAdataMultiBin <- redimMultiBinByClus(DNAdataMultiBinByClus)
    neighbourTree <- phangorn::rNNI(internalPhylo, moves = numMovesNNI)
    if (!identical(neighbourTree$tip.label, colnames(DNAdataMultiBin[[1]][[1]]))) {
        stop("Ordering of columns differ... Computing log-lik in internalPhylo... \n")
    } else{}
    updatedLogLik <- .logLikCpp(edgeMat = neighbourTree$edge, limProbsVec = limProbs, transMatList = transMatList, numOpenMP = numLikThreads, alignmentBin = DNAdataMultiBin, internalFlag = TRUE) ## Note that all edges belong to the same rate category, hence extTransMatList and intTransMatList taking the same value.
    MHratio <- exp(updatedLogLik - currentLogLik) ## Prior doesn't change...
    if (runif(1) < MHratio) {

        output <- list(internalPhylo = neighbourTree, logLik = updatedLogLik)
    } else{
        output <- list(internalPhylo = internalPhylo, logLik = currentLogLik)
    }
    output
}

.updateClusterPhylos <- function(currentValue, limProbs, intTransMatList, clusTransMatList, numMovesNNI = 1, numLikThreads, DNAdataBin, clusPhyloUpdateProp = 1) {

    logLikNow <- currentValue$logLik
    updateClusPhyloInt <- function(clusterPhylo) {

        if (!(class(clusterPhylo) == "phylo")) {
            output <- NA
        } else if (ape::Ntip(clusterPhylo) == 2) {
            output <- NA
        } else{
            output <- phangorn::rNNI(clusterPhylo, moves = numMovesNNI)
        }
        output
    }
    if (clusPhyloUpdateProp < 1) {
        treesSelectNum <- ceiling(length(currentValue$paraValues$clusterPhylos)*clusPhyloUpdateProp)
        selectIndices <- sample(seq_along(currentValue$paraValues$clusterPhylos), size = treesSelectNum, replace = FALSE) ## Uniformly selecting cluster-specific phylogenies for updates does not affect the transition kernel ratio. Also, I doubt it will greatly affect autocorrelation in the chain for cluster assignment indices.
    } else {
        selectIndices <- seq_along(currentValue$paraValues$clusterPhylos)
    }
    newTreesWithNAs <- lapply(currentValue$paraValues$clusterPhylos[selectIndices], FUN = updateClusPhyloInt)
    names(newTreesWithNAs) <- names(currentValue$paraValues$clusterPhylos[selectIndices])
    if (class(currentValue$paraValues$internalPhylo) != "phylo") { ## There's only one cluster... This is an ugly bit of code: improve ASAP.
      newLogLik <- .logLikCpp(edgeMat = newTreesWithNAs[[1]]$edge, limProbsVec = limProbs, transMatList = clusTransMatList, equivVector = names(limProbs), alignmentBin = lapply(DNAdataBin, FUN = function(x) {x[,currentValue$paraValues$clusterPhylos[[1]]$tip.label]}), returnRootMat = FALSE, internalFlag = FALSE, numOpenMP = numLikThreads)
      MHratio <- exp(newLogLik - logLikNow) ## Prior doesn't change...
      if (runif(1) < MHratio) {

        currentValue$logPostProb <- currentValue$logPostProb - currentValue$logLik + newLogLik
        currentValue$logLik <- newLogLik

        currentValue$DNAdataMultiBinByClus[[1]] <- .logLikCpp(edgeMat = newTreesWithNAs[[1]]$edge, limProbsVec = limProbs, transMatList = clusTransMatList, equivVector = names(limProbs), alignmentBin = DNAdataBin, returnRootMat = TRUE, internalFlag = FALSE, numOpenMP = numLikThreads)

        currentValue$paraValues$clusterPhylos[[1]] <- newTreesWithNAs[[1]] ## There's only one cluster, it follows that there are no NAs in newTreesWithNAs
        ## The prior doesn't change, so we don't need to update logPostProb.
      } else{}
      return(currentValue)
    } else{}
    NAbool <- sapply(newTreesWithNAs, FUN = function(x) class(x) != "phylo") ## NAs correspond to singletons and transmission pairs, since a NNI move does not change them.
    if (all(NAbool)) {
        return(currentValue) ## If no phylogenies can be updated, we simply return the current value of the chain.
    } else{}
    newTreesWithoutNAs <- newTreesWithNAs[!NAbool]
    alignmentBinList <- lapply(newTreesWithoutNAs, FUN = function(x) .dataBinSubset(DNAdataBin = DNAdataBin, keepSeqNamesOrNums = x$tip.label))
    edgeMatList <- lapply(newTreesWithoutNAs, FUN = function(x) x$edge)

    newClusMats <- .logLikCppV(edgeMatL = edgeMatList, limProbsVec = limProbs, transMatList = clusTransMatList, equivVector = names(limProbs), alignmentBinL = alignmentBinList, returnRootMat = TRUE, numOpenMP = numLikThreads, internalFlag = FALSE, priorBySizeTransMatBool = FALSE)
    names(newClusMats) <- names(newTreesWithoutNAs)

    internalClusPhylos <- function(x) {
        DNAdataMultiBinByClusCurr <- currentValue$DNAdataMultiBinByClus
        DNAdataMultiBinByClusCurr[[names(newClusMats)[x]]] <- newClusMats[[x]] ###################
        intAlignmentBin <- redimMultiBinByClus(DNAdataMultiBinByClusCurr)

        newLogLik <- .logLikCpp(edgeMat = currentValue$paraValues$internalPhylo$edge, limProbsVec = limProbs, transMatList = intTransMatList, equivVector = names(limProbs), alignmentBin = intAlignmentBin, returnRootMat = FALSE, internalFlag = TRUE, numOpenMP = numLikThreads)
        #MHratio <- exp(newLogLik - logLikNow) ## Prior doesn't change... Shouldn't it be currentValue$logLik
        MHratio <- exp(newLogLik - currentValue$logLik)
        if (runif(1) < MHratio) {

            labelClusToChange <- names(newClusMats)[x]
            currentValue$logPostProb <<- currentValue$logPostProb - currentValue$logLik + newLogLik
            currentValue$logLik <<- newLogLik

            currentValue$DNAdataMultiBinByClus <<- DNAdataMultiBinByClusCurr
            currentValue$paraValues$clusterPhylos[[labelClusToChange]] <<- newTreesWithoutNAs[[labelClusToChange]]
            ## The prior doesn't change, so we don't need to update logPostProb.
        } else{}
        NULL
    }
    lapply(seq_along(newClusMats), FUN = internalClusPhylos)

    currentValue
}

.logLikCppV <- function(edgeMatL, limProbsVec, alignmentMatL, transMatList, numOpenMP, equivVector, alignmentBinL, childNodeInClusIndic = 1000, returnRootMat = FALSE, internalFlag = FALSE, priorBySizeTransMatBool) {
    numStatesCons <- length(limProbsVec)
    ## The placeholders exist because Rcpp will want these arguments to be defined.
    placeholderAlignmentNum <- list(matrix(0,1,1))
    placeholderAlignmentAlpha <- list(matrix("a",1,1))
    placeholderBin <- list(list(matrix(1000,1,1)))
    placeholderEquiv <- rep("a",numStatesCons)
    if (missing(alignmentBinL)) {
        if (is.character(alignmentMatL[[1]][1,1])) {## Alphanumeric alignment, the function will need to convert it. It is fed to alignmentAlphaMat in logLikCppToWrap.
            output <- logLikCppToWrapV(edgeMatList = edgeMatL, alignmentMatList = placeholderAlignmentNum, limProbsVec = limProbsVec, transMatList = transMatList, numOpenMP = numOpenMP, equivVector = equivVector, alignmentAlphaMatList = alignmentMatL, alignmentBinList = placeholderBin, childNodeInClusIndic = childNodeInClusIndic, returnMatIndic = returnRootMat, internalFlag = internalFlag, priorBySizeTransMatBool = priorBySizeTransMatBool)
        } else {
            output <- logLikCppToWrapV(edgeMatList = edgeMatL, alignmentMatList = alignmentMatL, limProbsVec = limProbsVec, transMatList = transMatList, numOpenMP = numOpenMP, equivVector = placeholderEquiv, alignmentAlphaMatList = placeholderAlignmentAlpha, alignmentBinList = placeholderBin, childNodeInClusIndic = childNodeInClusIndic, returnMatIndic = returnRootMat, internalFlag = internalFlag, priorBySizeTransMatBool = priorBySizeTransMatBool) ## Both equivVector and alignmentAlphaMat are given arbitrary values serving as placeholders, since the underlying C++ function cannot handle missing arguments.
        }
    } else {
        output <- logLikCppToWrapV(edgeMatList = edgeMatL, alignmentMatList = placeholderAlignmentNum, limProbsVec = limProbsVec, transMatList = transMatList, numOpenMP = numOpenMP, equivVector = placeholderEquiv, alignmentAlphaMatList = placeholderAlignmentAlpha, alignmentBinList = alignmentBinL, childNodeInClusIndic = childNodeInClusIndic, returnMatIndic = returnRootMat, internalFlag = internalFlag, priorBySizeTransMatBool = priorBySizeTransMatBool)
    }
    output
}


.checkInput <- function(startingValues, Qmatrix, alignment, limProbs, shiftForAlpha) {
    ## Check if all necessary starting values are defined...
    test <- all(names(startingValues) %in% c("alpha", "clusInd", "phylogeny"))
    if (!test) {
        cat("Missing elements in startingValue list: ", paste(setdiff(c("alpha", "clusInd", "phylogeny"), names(startingValues)), sep = ", "),"\n")
        stop()
    } else{}
    ## Check if alignment has named rows...
    if (is.null(rownames(alignment))) { ## The tip ordering in the starting phylogeny should match the order of the rows in the alignment.
        warning("The alignment has unnamed rows (sequences). Make sure that the order of the tips in startingValues$phylogeny matches that of the rows in alignment. \n")
    } else{}

    ## Check if the alignment has as many tips as the starting phylogeny.
    test <- (nrow(alignment) == ape::Ntip(startingValues$phylogeny))
    if (!test) {
        stop("The number of rows (sequences) in the alignment is not equal to the number of tips in the starting value for the phylogeny. \n")
    } else{}
    ## Check if the names in clusInd are the same as in startingValues$phylogeny.
    test <- all(startingValues$phylogeny$tip.label %in% rownames(alignment)) & all(rownames(alignment) %in% startingValues$phylogeny$tip.label)
    if (!test) {
        stop("Mismatch between the row names in alignment and the tip labels in startingValues$phylogeny")
    } else{}
    ## Check if values in startingValues$clusInd correspond to clades in startingValues$phylogeny
    splitClusInd <- split(names(startingValues$clusInd), f = startingValues$clusInd)
    lapply(splitClusInd, FUN = function(x) {
        if (length(x) > 1) {
            mrcaNode <- ape::getMRCA(phy = startingValues$phylogeny, tip = x)
            childrenTips <- startingValues$phylogeny$tip.label[phangorn::Descendants(x = startingValues$phylogeny, node = mrcaNode, type = "tips")[[1]]]
            test <- setequal(childrenTips, x)
            if (!test) {
                stop("Clusters defined in startingValues$clusInd do not correspond to clades in startingValues$phylogeny.\n")
            } else{}
        } else{}
        NULL
    })

    ## Checking if limProbs has named elements. It is essential for handling ambiguities!
    test <- !is.null(names(limProbs))
    if (!test) {
        stop("The limProbs vector must have named elements. (Also make sure that its ordering matches that of rows in Qmatrix.)\n")
    } else{}

    ## Qmatrix should have named rows...
    if (!missing(Qmatrix)) {
      test <- !is.null(rownames(Qmatrix))
      if (!test) {
        warning("The Q matrix does not have row names. Assuming that the row order matches that of the elements of limProbs... \n")
      } else {
        newTest <- setequal(rownames(Qmatrix), names(limProbs))
        if (!newTest) {
          stop("The names of the elements in limProbs do not match the row names of Qmatrix.\n")
        } else{}
      }
    } else{}
    ## There should not be any capital letters in a nucleotide alignment: it prevents identification of ambiguities.
    test <- !any(c("A", "T", "C", "G") %in% alignment)
    if (!test) {
        warning("alignment contains capital letters. Is the alignment made up of nucleotides? If so, be advised that using small letters - i.e. a,c,t, and g - lets the algorithm recognize and correctly handle ambiguities. \n")
    } else{}
    ## The shift parameter for alpha should be positive.
    test <- shiftForAlpha >= 0
    if (!test) {
        stop("shiftForAlpha should take value 0 or greater.\n")
    } else{}
    NULL
}

.checkArgumentsLogLikFromSplitPhylo <- function() {}

.checkArgumentsLogLikFromClusInd <- function() {}



