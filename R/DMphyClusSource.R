.performStepPhylo <- function(currentValue, limProbs, shapePriorAlpha, scalePriorAlpha, numGammaCat, betweenTransMatAll, withinTransMatAll, enableJainNeal = FALSE, currentIter, DNAdata, numNeighbours, numMovesNNIbetween = 1, numMovesNNIwithin = 1, numLikThreads, DNAdataBin, poisRateNumClus, clusPhyloUpdateProp = 1, alphaMin, numSplitMergeMoves) {
    if (max(currentValue$paraValues$clusInd) > 1) {
      currentValue <- .updateTransMatrix(currentValue = currentValue, allTransMatList = betweenTransMatAll, transMatsNoChange = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], samplingRadius = 2, betweenBool = TRUE, limProbs = limProbs, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin)
    }
    currentValue <- .updateTransMatrix(currentValue = currentValue, allTransMatList = withinTransMatAll, transMatsNoChange = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], samplingRadius = 2, betweenBool = FALSE, limProbs = limProbs, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin)

    if (length(currentValue$paraValues$clusterPhylos) > 2) {
        internalPhyloList <- .updateSupportingPhylo(internalPhylo = currentValue$paraValues$internalPhylo, logLimProbs = logLimProbs, logTransMatList = logIntTransMatAll[[currentValue$paraValues$logIntMatListIndex]], numMovesNNI = numMovesNNIint, numLikThreads = numLikThreads, DNAdataMultiBinByClus = currentValue$DNAdataMultiBinByClus, currentLogLik = currentValue$logLik) ## This step updates the supporting phylogeny.
        currentValue$logPostProb <- currentValue$logPostProb - currentValue$logLik + internalPhyloList$logLik ## Updating the internal phylogeny doesn't change the prior value.
        currentValue$logLik <- internalPhyloList$logLik
        currentValue$paraValues$internalPhylo <- internalPhyloList$internalPhylo
    } else{}

    currentValue <- .updateClusterPhylos(currentValue = currentValue, logLimProbs = logLimProbs, intTransMatList = logIntTransMatAll[[currentValue$paraValues$logIntMatListIndex]], clusTransMatList = logExtTransMatAll[[currentValue$paraValues$logExtMatListIndex]], numMovesNNI = numMovesNNIext, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin, clusPhyloUpdateProp = clusPhyloUpdateProp) ## This step updates the cluster-specific phylogenies. It is the slowest step.
    lapply(1:numSplitMergeMoves, FUN = function(x) { ##
        currentValue <<- .splitJoinClusterMove(currentValue = currentValue, DNAdataBin = DNAdataBin, logLimProbs = logLimProbs, clusTransMatList = logExtTransMatAll[[currentValue$paraValues$logExtMatListIndex]], intTransMatList = logIntTransMatAll[[currentValue$paraValues$logIntMatListIndex]], numLikThreads = numLikThreads, singletonMatrices = singletonMatrices, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapePriorAlpha, scaleForAlpha = scalePriorAlpha, alphaMin = alphaMin)
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

.initCurrentValue <- function(startingValues, DNAdataBin, betweenTransMatAll, withinTransMatAll, limProbs, numLikThreads, shapeForAlpha, scaleForAlpha, alphaMin) {
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
    currentValue$paraValues$betweenMatListIndex <- ceiling(length(betweenTransMatAll)/2)
    currentValue$paraValues$withinMatListIndex <- ceiling(length(withinTransMatAll)/2)
    currentValue$paraValues$clusInd <- startingValues$clusInd[startingValues$phylogeny$tip.label]

    if (max(currentValue$paraValues$clusInd) == 1) {
      currentValue$paraValues$clusterNodeIndices <- ape::Ntip(startingValues$phylogeny)
      names(currentValue$paraValues$clusterNodeIndices) <- "1"
    } else {
      currentValue$paraValues$clusterNodeIndices <- sapply(1:max(currentValue$paraValues$clusInd), FUN = function(clusterIndex) {
        seqsInCluster <- names(currentValue$paraValues$clusInd)[currentValue$paraValues$clusInd == clusterIndex]
        if (length(seqsInCluster) > 1) {
          return(ape::getMRCA(currentValue$paraValues$phylogeny, seqsInCluster))
        }
        match(seqsInCluster, currentValue$paraValues$tip.label)
      })
    }
    currentValue$logLik <- logLikCpp(edgeMat = currentValue$paraValues$phylogeny$edge, limProbsVec = limProbs, withinTransMatList = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], betweenTransMatList = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], numOpenMP = numLikThreads, alignmentBin = DNAdataBin, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, numTips = ape::Ntip(currentValue$paraValues$phylogeny))$logLik
    currentValue$clusterCounts <- as.vector(table(currentValue$paraValues$clusInd)) ## The as.vector ensures that clusterCounts behaves always as a vector, but it removes the names.
    names(currentValue$clusterCounts) <- 1:max(currentValue$paraValues$clusInd) ## This once again rests on the assumption that there is no gap in the cluster labels. This is an assumption we made before. Names are needed sometimes, although it is better to index by number, rather than by name (must require looking into an index).
    currentValue$logPostProb <- currentValue$logLik + clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE) ## The starting value for log-posterior probability for the starting partition.
    currentValue
}

.performSplit <- function(currentValue, numSplits, clusLabel, DNAdataBin, logLimProbs, clusTransMatList, intTransMatList, numLikThreads, singletonMatrices, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    originalPhylo <- currentValue$paraValues$clusterPhylos[[clusLabel]]
    newClusLabel <- as.character(max(currentValue$paraValues$clusInd)+1)
    graftedTree <- ape::rtree(2)
    graftedTree$edge.length <- NULL
    graftedTree$tip.label <- c(clusLabel, newClusLabel)
    if (length(currentValue$paraValues$clusterPhylos) > 1) {
        newInternalTree <- suppressWarnings(ape::bind.tree(x = ape::reorder.phylo(currentValue$paraValues$internalPhylo), y = graftedTree, where = match(clusLabel, currentValue$paraValues$internalPhylo$tip.label))) ## Breaking up a cluster results in the internal phylogeny changing too, since it must now integrate the part of the cluster phylogeny that used to link the two new clusters. The ordering of the labels in the size 2 phylogeny supporting the clusters doesn't matter: it is symmetrical.
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
            newSitePatterns <- getSitePatterns(alignmentBin)
            if (!identical(colnames(alignmentBin[[1]]), newSubTree$tip.label)) {
                stop("Look at names... \n")
            } else{}
            return(list(tree = newSubTree, clusMatrix = .logLikCpp(edgeMat = newSubTree$edge, logLimProbsVec = logLimProbs, logTransMatList = clusTransMatList, numOpenMP = numLikThreads, alignmentBin = newSitePatterns$uniqueDNAdataBin, internalFlag = FALSE, returnRootMat = TRUE, sitePatterns = newSitePatterns$sitePatterns), sitePatterns = newSitePatterns))
        }
    })
    DNAdataMultiBinByClusCopy[graftedTree$tip.label] <- lapply(splitDNAbyClus, FUN = function(x) x$clusMatrix)
    newClusterPhylo <- lapply(splitDNAbyClus, FUN = function(x) x$tree)
    newSitePatterns <- lapply(splitDNAbyClus, FUN = function(x) x$sitePatterns)
    newClusPhylosLabel <- graftedTree$tip.label

    DNAdataMultiBinByClusCopy <- DNAdataMultiBinByClusCopy[newInternalTree$tip.label]

    if (!identical(names(DNAdataMultiBinByClusCopy), newInternalTree$tip.label)) {
        stop("Look at element names! \n")
    } else{}
    newLogLik <- .logLikCpp(edgeMat = newInternalTree$edge, logLimProbsVec = logLimProbs, logTransMatList = intTransMatList, numOpenMP = numLikThreads, alignmentBin = DNAdataMultiBinByClusCopy, internalFlag = TRUE, returnRootMat = FALSE)
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
    ancestorsGroupsPairsNew <- ancestorsGroupsNew[vapply(ancestorsGroupsNew, FUN = function(x) {length(x) > 1}, FUN.VALUE = c(Return = TRUE))]
    numPairsNew <- length(ancestorsGroupsPairsNew)
    transKernRatio <- numSplits/numPairsNew
    list(logLik = newLogLik, clusInd = newClusInd, counts = newCounts, logPostProb = newLogPostProb, transKernRatio = transKernRatio, clusterPhylos = newClusterPhylo, internalPhylo = newInternalTree, DNAdataMultiBinByClus = DNAdataMultiBinByClusCopy, newClusPhylosLabel = newClusPhylosLabel, sitePatterns = newSitePatterns)
}

.performMerge <- function(currentValue, ancestorsGroupsPairs, pairToSelect, DNAdataBin, logLimProbs, clusTransMatList, intTransMatList, numLikThreads, tipAncestors, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    numPairs <- length(ancestorsGroupsPairs)
    clustersInPair <- names(ancestorsGroupsPairs[[pairToSelect]])
    hostTree <- ape::rtree(2, tip.label = c(0,0))
    hostTree$edge.length <- 0
    if (class(currentValue$paraValues$clusterPhylos[[clustersInPair[[1]]]]) == "phylo") {
        hostTree <- suppressWarnings(ape::bind.tree(x = hostTree, y = ape::reorder.phylo(currentValue$paraValues$clusterPhylos[[clustersInPair[[1]]]]), where = 1))
        hostChildren <- phangorn::Children(hostTree, ape::Ntip(hostTree)+1)
        whereToBind <- hostChildren[hostChildren < ape::Ntip(hostTree)]
    } else{
        hostTree$tip.label[1] <- names(currentValue$paraValues$clusInd)[match(as.numeric(clustersInPair[[1]]), currentValue$paraValues$clusInd)]
        whereToBind <- 2
    }

    if (class(currentValue$paraValues$clusterPhylos[[clustersInPair[[2]]]]) == "phylo") {
        newClusterPhylo <- suppressWarnings(ape::bind.tree(x = hostTree, y = ape::reorder.phylo(currentValue$paraValues$clusterPhylos[[clustersInPair[[2]]]]), where = whereToBind))
    } else{
        hostTree$tip.label[whereToBind] <- names(currentValue$paraValues$clusInd)[match(as.numeric(clustersInPair[[2]]), currentValue$paraValues$clusInd)]
        newClusterPhylo <- hostTree
    }
    newTipLabel <- as.character(max(currentValue$paraValues$clusInd)+1)
    currentValue$paraValues$internalPhylo$node.label[[tipAncestors[[clustersInPair[[1]]]]-ape::Ntip(currentValue$paraValues$internalPhylo)]] <- newTipLabel ## The new cluster will be given the label of a totally new cluster.
    newInternalTree <- suppressWarnings(ape::drop.tip(currentValue$paraValues$internalPhylo, tip = clustersInPair, trim.internal = FALSE)) ## The node label should become a tip label. There will be a gap in the nodes. We suppress warnings, because the algorithm can handle the case where all tips are dropped.
    mergeDataBin <- .dataBinSubset(DNAdataBin = DNAdataBin, keepSeqNamesOrNums = newClusterPhylo$tip.label)
    newSitePatterns <- getSitePatterns(mergeDataBin)
    DNAdataMultiBinByClusCopy <- currentValue$DNAdataMultiBinByClus[-match(clustersInPair, names(currentValue$DNAdataMultiBinByClus))]
    clusTransMatSelectIndex <- min(length(clusTransMatList), ape::Ntip(newClusterPhylo))
    if (!identical(newClusterPhylo$tip.label, colnames(mergeDataBin[[1]]))) {
        stop("Look at colnames! \n")
    } else{}
    DNAdataMultiBinByClusCopy[[newTipLabel]] <- .logLikCpp(edgeMat = newClusterPhylo$edge, logLimProbsVec = logLimProbs, logTransMatList = clusTransMatList, numOpenMP = numLikThreads, alignmentBin = newSitePatterns$uniqueDNAdataBin, internalFlag = FALSE, returnRootMat = TRUE, sitePatterns = newSitePatterns$sitePatterns) ## This is automatically added to the end of the list, which means the ordering of the list should not be affected.

    newClusPhylosLabel <- newTipLabel

    if (!is.null(newInternalTree)) {
      if (!identical(newInternalTree$tip.label, names(DNAdataMultiBinByClusCopy))) {
        stop("Look at ordering of DNAdataMultiBinByClusCopy! \n")
      } else{}
      newLogLik <- .logLikCpp(edgeMat = newInternalTree$edge, logLimProbsVec = logLimProbs, logTransMatList = intTransMatList, numOpenMP = numLikThreads, alignmentBin = DNAdataMultiBinByClusCopy, internalFlag = TRUE, returnRootMat = FALSE)
    } else { ## If the internal tree is NULL, it means we only have one cluster!
      if (!identical(colnames(mergeDataBin[[1]]), newClusterPhylo$tip.label)) {
        stop("Error in .performMerge: orders of columns in mergeDataBin should match that of the tip labels in newClusterPhylo!\n")
      } else{}
      newLogLik <- .logLikCpp(edgeMat = newClusterPhylo$edge, logLimProbsVec = logLimProbs, logTransMatList = clusTransMatList, numOpenMP = numLikThreads, alignmentBin = newSitePatterns$uniqueDNAdataBin, internalFlag = FALSE, returnRootMat = FALSE, sitePatterns = newSitePatterns$sitePatterns) ## Should make sure that the order in mergeDataBin matches that in newClusterPhylo...
    }
    newClusInd <- replace(currentValue$paraValues$clusInd, which(currentValue$paraValues$clusInd %in% as.numeric(clustersInPair)), as.numeric(newTipLabel)) ## There's a gap in the cluster labels.
    newCounts <- table(newClusInd)
    newLogPostProb <- newLogLik + clusIndLogPrior(clusInd = newClusInd, alpha = currentValue$paraValues$alpha) + dpois(length(newCounts), lambda = poisRateNumClus, log = TRUE) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)## Added the Poisson log-prob. to reflect a Poisson prior on the total number of clusters.
    transKernRatio <- numPairs/sum(newCounts>1)
    newClusterPhylo <- list(newClusterPhylo)
    newSitePatterns <- list(newSitePatterns)
    list(logLik = newLogLik, clusInd = newClusInd, counts = newCounts, logPostProb = newLogPostProb, transKernRatio = transKernRatio, clusterPhylos = newClusterPhylo, internalPhylo = newInternalTree, DNAdataMultiBinByClus = DNAdataMultiBinByClusCopy, newClusPhylosLabel = newClusPhylosLabel, sitePatterns = newSitePatterns)
}

.splitJoinClusterMove <- function(currentValue, DNAdataBin, logLimProbs, clusTransMatList, intTransMatList, numLikThreads, singletonMatrices, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    if (class(currentValue$paraValues$internalPhylo) == "phylo") {
      tipAncestors <- phangorn::Ancestors(currentValue$paraValues$internalPhylo, node = seq_along(currentValue$paraValues$internalPhylo$tip.label), type = "parent")
      names(tipAncestors) <- currentValue$paraValues$internalPhylo$tip.label
      ancestorsGroups <- split(tipAncestors, f = tipAncestors)
      ancestorsGroupsPairs <- ancestorsGroups[vapply(ancestorsGroups, FUN = function(x) {length(x) > 1}, FUN.VALUE = logical(1))]
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
        splitMergeResult <- .performSplit(currentValue = currentValue, numSplits = length(clustersToSplit), clusLabel = clusLabel, DNAdataBin = DNAdataBin, logLimProbs = logLimProbs, clusTransMatList = clusTransMatList, intTransMatList = intTransMatList, numLikThreads = numLikThreads, singletonMatrices = singletonMatrices, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
    } else {
        pairToSelect <- sample.int(n = numPairs, size = 1)
        clustersInPair <- names(ancestorsGroupsPairs[[pairToSelect]])
        splitMergeResult <- .performMerge(currentValue = currentValue, ancestorsGroupsPairs = ancestorsGroupsPairs, pairToSelect = pairToSelect, DNAdataBin = DNAdataBin, logLimProbs = logLimProbs, clusTransMatList = clusTransMatList, intTransMatList = intTransMatList, numLikThreads = numLikThreads, tipAncestors = tipAncestors, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
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
        currentValue$sitePatternsByClus[splitMergeResult$newClusPhylosLabel] <- splitMergeResult$sitePatterns

        if (!splitMove) {

            matchingClusters <- match(clustersInPair, names(currentValue$paraValues$clusterPhylos))
            currentValue$paraValues$clusterPhylos <- currentValue$paraValues$clusterPhylos[-matchingClusters] ## Removing the clusters that have been merged into a new cluster.
            names(currentValue$paraValues$clusterPhylos) <- as.character(.relabel(as.numeric(names(currentValue$paraValues$clusterPhylos))))
            currentValue$sitePatternsByClus <- currentValue$sitePatternsByClus[-matchingClusters]
            names(currentValue$sitePatternsByClus) <- as.character(.relabel(as.numeric(names(currentValue$sitePatternsByClus))))
            names(currentValue$DNAdataMultiBinByClus) <- as.character(.relabel(as.numeric(names(currentValue$DNAdataMultiBinByClus))))
            names(currentValue$clusterCounts) <- as.character(.relabel(as.numeric(names(currentValue$clusterCounts))))
            currentValue$paraValues$clusInd <- .relabel(currentValue$paraValues$clusInd)
            currentValue$paraValues$internalPhylo$tip.label <- as.character(.relabel(as.numeric(currentValue$paraValues$internalPhylo$tip.label)))
        } else {}
    } else{}
    currentValue
}

.updateTransMatrix <- function(currentValue, allTransMatList, transMatsNoChange, samplingRadius = 2, betweenBool, limProbs, numLikThreads, DNAdataBin) {
  if (betweenBool) {
    basicIndex <- currentValue$paraValues$betweenMatListIndex
  } else{
    basicIndex <- currentValue$paraValues$withinMatListIndex
  }
  newIndex <- basicIndex + sample(c(-samplingRadius:-1,1:samplingRadius), size = 1)
  newIndex <- (newIndex<=0)*(length(allTransMatList) + newIndex) + (newIndex > length(allTransMatList))*(newIndex - length(allTransMatList)) + ((newIndex > 0) & (newIndex <= length(allTransMatList)))*newIndex  ## We have a circular transition kernel. If the range for allTransMatList is wide enough, the circularity should be never invoked.
  if (betweenBool) {
    betweenTransMatList <- allTransMatList[[newIndex]]
    withinTransMatList <- transMatsNoChange
  } else {
    betweenTransMatList <- transMatsNoChange
    withinTransMatList <- allTransMatList[[newIndex]]
  }
  newLogLik <- logLikCpp(edgeMat = currentValue$paraValues$phylogeny$edge, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, limProbsVec = limProbs, withinTransMatList = withinTransMatList, betweenTransMatList = betweenTransMatList, numOpenMP = numLikThreads, alignmentBin = DNAdataBin, numTips = ape::Ntip(currentValue$paraValues$phylogeny))
 
  MHratio <- exp(newLogLik$logLik - currentValue$logLik)
  if (runif(1) < MHratio) {
    if (betweenBool) {
      currentValue$paraValues$betweenMatListIndex <- newIndex
    } else{
      currentValue$paraValues$withinMatListIndex <- newIndex
    }
    currentValue$logPostProb <- newLogLik$logLik +  (currentValue$logPostProb - currentValue$logLik)
    currentValue$logLik <- newLogLik$logLik
  } else{}
  currentValue
}

.DMphyClusCore <- function(nIter, startingValues, limProbs, shapeForAlpha, scaleForAlpha, numMovesNNIbetween, numMovesNNIwithin, numLikThreads, DNAdataBin, poisRateNumClus, clusPhyloUpdateProp, numSplitMergeMoves, alphaMin, withinTransMatAll, betweenTransMatAll) {

  if (is.matrix(withinTransMatAll[[1]])) {
    numGammaCat <- length(withinTransMatAll) ## We only have one set of substitution rate matrices.
  } else{
    numGammaCat <- length(withinTransMatAll[[1]])
  }

  currentValue <- .initCurrentValue(startingValues = startingValues, DNAdataBin = DNAdataBin, withinTransMatAll = withinTransMatAll, betweenTransMatAll = betweenTransMatAll, limProbs = limProbs, numLikThreads = numLikThreads, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)

  cat("Launching chain... \n \n")

  progress <- txtProgressBar(style = 3)
  setTxtProgressBar(pb = progress, value = 0)
  onePerc <- nIter/100
  longOut <- lapply(1:nIter, FUN = function(x) {
    if ((x %% onePerc) == 0) {
      setTxtProgressBar(pb = progress, value = x/nIter)
    } else{}

    currentValue <<- .performStepPhylo(currentValue = currentValue, limProbs = limProbs, shapePriorAlpha = shapeForAlpha, scalePriorAlpha = scaleForAlpha, withinTransMatAll = withinTransMatAll, betweenTransMatAll = betweenTransMatAll, currentIter = x, numMovesNNIbetween = numMovesNNIbetween, numMovesNNIwithin = numMovesNNIwithin, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin, poisRateNumClus = poisRateNumClus, clusPhyloUpdateProp = clusPhyloUpdateProp, alphaMin = alphaMin, numSplitMergeMoves = numSplitMergeMoves)

    paraVec <- currentValue$paraValues
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

.logLikCpp <- function(edgeMat, logLimProbsVec, logTransMatList, numOpenMP, equivVector, alignmentBin, returnRootMat = FALSE, internalFlag = FALSE, sitePatterns) {

    if (missing(sitePatterns) & !internalFlag) {
      sitePatterns <- seq_along(alignmentBin)
    } else if (internalFlag) {
      sitePatterns <- 1
    } else{}
    numStatesCons <- length(logLimProbsVec)
    placeholderEquiv <- rep("a",numStatesCons)
    output <- logLikCppToWrap(edgeMat = edgeMat, logLimProbsVec = logLimProbsVec, logTransMatList = logTransMatList, numOpenMP = numOpenMP, equivVector = placeholderEquiv, alignmentBin = alignmentBin, returnMatIndic = returnRootMat, internalFlag = internalFlag, sitePatterns = sitePatterns)

    output
}

.introduceMultiPhyloWithDist <-  function(phylogeny, clusInd) {
    if (is.null(phylogeny$tip.label)) {
        toSplit <- as.character(1:ape::Ntip(phylogeny))
    } else {
        toSplit <- phylogeny$tip.label
    }
    nodeListCut <- split(x = toSplit, f = clusInd)
    keepIndices <- vapply(nodeListCut, FUN = function(x) {length(x) > 2}, FUN.VALUE = c(Logical = TRUE))
    nodeListCut <- nodeListCut[keepIndices]
    newPhylo <- phylogeny
    internalFun <- function(seqsInClusInd) {
        multiTips <- which(newPhylo$tip.label %in% seqsInClusInd)
        mrcaNode <- ape::getMRCA(phy = newPhylo, tip = multiTips)
        distsToMRCA <- as.vector(ape::dist.nodes(newPhylo)[mrcaNode, match(seqsInClusInd, newPhylo$tip.label)])
        newPhylo$tip.label[multiTips] <- "removeMe"
        newTree <- ape::stree(length(seqsInClusInd), tip.label = seqsInClusInd)
        newTree$edge.length <- distsToMRCA
        incrementedTree <- ape::reorder.phylo(suppressWarnings(ape::bind.tree(x = newPhylo, y = newTree, where = mrcaNode)))
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
.updateSupportingPhylo <- function(internalPhylo, currentLogLik, logLimProbs, logTransMatList, numMovesNNI, numLikThreads, DNAdataMultiBinByClus) {

    neighbourTree <- phangorn::rNNI(internalPhylo, moves = numMovesNNI)
    if (!identical(neighbourTree$tip.label, names(DNAdataMultiBinByClus))) {
        stop("Ordering of columns differ... Computing log-lik in internalPhylo... \n")
    } else{}
    updatedLogLik <- .logLikCpp(edgeMat = neighbourTree$edge, logLimProbsVec = logLimProbs, logTransMatList = logTransMatList, numOpenMP = numLikThreads, alignmentBin = DNAdataMultiBinByClus, internalFlag = TRUE) ## Note that all edges belong to the same rate category, hence extTransMatList and intTransMatList taking the same value.
    MHratio <- exp(updatedLogLik - currentLogLik) ## Prior doesn't change...
    if (runif(1) < MHratio) {

        output <- list(internalPhylo = neighbourTree, logLik = updatedLogLik)
    } else{
        output <- list(internalPhylo = internalPhylo, logLik = currentLogLik)
    }
    output
}

.updateClusterPhylos <- function(currentValue, logLimProbs, intTransMatList, clusTransMatList, numMovesNNI = 1, numLikThreads, DNAdataBin, clusPhyloUpdateProp = 1) {

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
      newLogLik <- .logLikCpp(edgeMat = newTreesWithNAs[[1]]$edge, logLimProbsVec = logLimProbs, logTransMatList = clusTransMatList, equivVector = names(logLimProbs), alignmentBin = currentValue$sitePatternsByClus[[1]]$uniqueDNAdataBin, returnRootMat = FALSE, internalFlag = FALSE, numOpenMP = numLikThreads, sitePatterns = currentValue$sitePatternsByClus[[1]]$sitePatterns)
      MHratio <- exp(newLogLik - logLikNow) ## Prior doesn't change...
      if (runif(1) < MHratio) {

        currentValue$logPostProb <- currentValue$logPostProb - currentValue$logLik + newLogLik
        currentValue$logLik <- newLogLik

        currentValue$DNAdataMultiBinByClus[[1]] <- .logLikCpp(edgeMat = newTreesWithNAs[[1]]$edge, logLimProbsVec = logLimProbs, logTransMatList = clusTransMatList, equivVector = names(logLimProbs), alignmentBin = currentValue$sitePatternsByClus[[1]]$uniqueDNAdataBin, returnRootMat = TRUE, internalFlag = FALSE, numOpenMP = numLikThreads, sitePatterns = currentValue$sitePatternsByClus[[1]]$sitePatterns)

        currentValue$paraValues$clusterPhylos[[1]] <- newTreesWithNAs[[1]] ## There's only one cluster, it follows that there are no NAs in newTreesWithNAs
        ## The prior doesn't change, so we don't need to update logPostProb.
      } else{}
      return(currentValue)
    } else{}
    NAbool <- vapply(newTreesWithNAs, FUN = function(x) {class(x) != "phylo"}, FUN.VALUE = logical(1)) ## NAs correspond to singletons and transmission pairs, since a NNI move does not change them.
    if (all(NAbool)) {
        return(currentValue) ## If no phylogenies can be updated, we simply return the current value of the chain.
    } else{}
    newTreesWithoutNAs <- newTreesWithNAs[!NAbool]

    newClusMats <- lapply(names(newTreesWithoutNAs), FUN = function(phyloName) {
      edgeMat <- newTreesWithoutNAs[[phyloName]]$edge
      alignmentBin <- currentValue$sitePatternsByClus[[phyloName]]$uniqueDNAdataBin
      sitePatterns <- currentValue$sitePatternsByClus[[phyloName]]$sitePatterns
      .logLikCpp(edgeMat = edgeMat, logLimProbsVec = logLimProbs, logTransMatList = clusTransMatList, equivVector = names(logLimProbs), alignmentBin = alignmentBin, returnRootMat = TRUE, numOpenMP = numLikThreads, internalFlag = FALSE, sitePatterns = sitePatterns)
    })
    names(newClusMats) <- names(newTreesWithoutNAs)

    internalClusPhylos <- function(x) {

        DNAdataMultiBinByClusCurr <- currentValue$DNAdataMultiBinByClus
        DNAdataMultiBinByClusCurr[[names(newClusMats)[x]]] <- newClusMats[[x]]
        newLogLik <- .logLikCpp(edgeMat = currentValue$paraValues$internalPhylo$edge, logLimProbsVec = logLimProbs, logTransMatList = intTransMatList, equivVector = names(logLimProbs), alignmentBin = DNAdataMultiBinByClusCurr, returnRootMat = FALSE, internalFlag = TRUE, numOpenMP = numLikThreads)
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

# getSitePatterns <- function(alignment) { ## Works for DNAdataBin type alignments by also for conventional alignments in matrix form (1 row for each sequence, 1 col for each locus)
#   if (is.list(alignment)) {
#     hashesVec <- vapply(alignment, FUN = digest::digest, algo = "xxhash64", FUN.VALUE = character(1))
#   } else {
#     hashesVec <- vapply(1:ncol(alignment), FUN = function(x) {digest::digest(alignment[,x], algo = "xxhash64")}, FUN.VALUE = character(1))
#   }
#   noDuplicates <- !duplicated(hashesVec)
#   uniqueHashes <- hashesVec[noDuplicates]
#   matchVector <- match(hashesVec, uniqueHashes)
#
#   if (is.list(alignment)) {
#     uniqueMat <- alignment[noDuplicates]
#   } else{
#     uniqueMat <- alignment[,noDuplicates]
#   }
#
#   list(uniqueDNAdataBin = uniqueMat, sitePatterns = matchVector)
# }
