.performStepPhylo <- function(currentValue, limProbs, shapePriorAlpha, scalePriorAlpha, numGammaCat, betweenTransMatAll, withinTransMatAll, enableJainNeal = FALSE, currentIter, DNAdata, numNeighbours, numMovesNNIbetween = 1, numMovesNNIwithin = 1, numLikThreads, DNAdataBin, poisRateNumClus, clusPhyloUpdateProp = 1, alphaMin, numSplitMergeMoves, maxMapSize = NULL, cullProportion = NULL, togglePhyloUpdates = TRUE) {
  
  if (max(currentValue$paraValues$clusInd) > 1) {
      currentValue <- .updateTransMatrix(currentValue = currentValue, allTransMatList = betweenTransMatAll, transMatsNoChange = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], samplingRadius = 2, betweenBool = TRUE, limProbs = limProbs, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin)
    }
    
    currentValue <- .updateTransMatrix(currentValue = currentValue, allTransMatList = withinTransMatAll, transMatsNoChange = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], samplingRadius = 2, betweenBool = FALSE, limProbs = limProbs, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin)
    if (togglePhyloUpdates)
    {
      if (max(currentValue$paraValues$clusInd) > 2) {
        currentValue <- .updateBetweenPhylo(currentValue = currentValue, limProbs = limProbs, withinTransMatList = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], betweenTransMatList = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], numMovesNNI = numMovesNNIbetween, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin) ## This step updates the supporting phylogeny.
      } else{}
      
      currentValue <- .updateClusterPhylos(currentValue = currentValue, limProbs = limProbs, betweenTransMatList = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], withinTransMatList = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], numMovesNNI = numMovesNNIwithin, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin, clusPhyloUpdateProp = clusPhyloUpdateProp)    ## This step updates the cluster-specific phylogenies. It is the slowest step.
    }
    lapply(1:numSplitMergeMoves, FUN = function(x) { ## Should this be made into a recursive function call?
        currentValue <<- .splitJoinClusterMove(currentValue = currentValue, DNAdataBin = DNAdataBin, limProbs = limProbs, withinTransMatList = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], betweenTransMatList = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], numLikThreads = numLikThreads, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapePriorAlpha, scaleForAlpha = scalePriorAlpha, alphaMin = alphaMin)
    })

    currentValue <- .updateAlpha(currentValue, shapePriorAlpha, scalePriorAlpha, alphaMin = alphaMin)
    if (!is.null(maxMapSize))
    {
      checkAndCullMap(AugTreePointer = currentValue$extPointer, allowedNumberOfElements = maxMapSize, cullProportion = cullProportion, withinTransProbs = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], betweenTransProbs = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], limProbs = limProbs)
    }
    currentValue
}

reorderTips <- function(phylogeny, newTipOrder)
{
  reordering <- match(phylogeny$tip.label, newTipOrder)
  equivalenceVec <- 1:max(phylogeny$edge)
  equivalenceVec[1:ape::Ntip(phylogeny)] <- reordering
  phylogeny$edge[,2] <- equivalenceVec[phylogeny$edge[,2]]
  phylogeny$tip.label <- newTipOrder
  phylogeny
}

.updateAlpha <- function(currentValue, shapeForAlpha, scaleForAlpha, alphaMin = 6) {
    increment <- (runif(1) - 0.5)
    newAlpha <- currentValue$paraValues$alpha + increment
    newAlpha <- newAlpha + abs(newAlpha - alphaMin)*(newAlpha < alphaMin)
    newLogPrior <- clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = newAlpha) + dgamma(newAlpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)
    currentLogPrior <- clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE)
    poisContrib <- currentValue$logPostProb - currentLogPrior
    
    MHratio <- exp(newLogPrior - currentLogPrior)
    if (runif(1) < MHratio) {
        currentValue$paraValues$alpha <- newAlpha
        currentValue$logPostProb <- currentValue$logLik + newLogPrior + poisContrib
    } else{}
    currentValue
}

.initCurrentValue <- function(startingValues, DNAdataBin, betweenTransMatAll, withinTransMatAll, limProbs, numLikThreads, shapeForAlpha, scaleForAlpha, alphaMin, poisRateNumClus) {
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
    currentValue$paraValues$betweenMatListIndex <- startingValues$betweenMatListIndex
    if (is.null(startingValues$betweenMatListIndex))
    {
      currentValue$paraValues$betweenMatListIndex <- ceiling(length(betweenTransMatAll)/2)
    }
    currentValue$paraValues$withinMatListIndex <- startingValues$withinMatListIndex
    if (is.null(startingValues$withinMatListIndex))
    {
      currentValue$paraValues$withinMatListIndex <- ceiling(length(withinTransMatAll)/2)
    }
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
        match(seqsInCluster, currentValue$paraValues$phylogeny$tip.label)
      })
    }
    logLikAndPointer <- logLikCpp(edgeMat = currentValue$paraValues$phylogeny$edge, limProbsVec = limProbs, withinTransMatList = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], betweenTransMatList = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], numOpenMP = numLikThreads, alignmentBin = DNAdataBin, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, numTips = ape::Ntip(currentValue$paraValues$phylogeny), numLoci = length(DNAdataBin[[1]]), withinMatListIndex = currentValue$paraValues$withinMatListIndex, betweenMatListIndex = currentValue$paraValues$betweenMatListIndex)
    currentValue$logLik <- logLikAndPointer$logLik
    currentValue$extPointer <- logLikAndPointer$solutionPointer
    currentValue$alignmentBinPointer <- logLikAndPointer$alignmentBinPointer
    currentValue$clusterCounts <- as.vector(table(currentValue$paraValues$clusInd)) ## The as.vector ensures that clusterCounts behaves always as a vector, but it removes the names.
    names(currentValue$clusterCounts) <- 1:max(currentValue$paraValues$clusInd) ## This once again rests on the assumption that there is no gap in the cluster labels. This is an assumption we made before. Names are needed sometimes, although it is better to index by number, rather than by name (must require looking into an index).
    currentValue$logPostProb <- currentValue$logLik + clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE) + dpois(max(currentValue$paraValues$clusInd), lambda = poisRateNumClus, log = TRUE) ## The starting value for log-posterior probability for the starting partition.
    currentValue
}

.performSplit <- function(currentValue, numMoves, clusNumber, DNAdataBin, limProbs, withinTransMatList, betweenTransMatList, numLikThreads, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

  currentPhylo <- currentValue$paraValues$phylogeny ## This is merely an alias. It won't be modified, so it won't be copied.
  
  newClusNodes <- phangorn::Children(currentPhylo, currentValue$paraValues$clusterNodeIndices[[clusNumber]]) 
  newClusMRCAs <- currentValue$paraValues$clusterNodeIndices
  newClusMRCAs[clusNumber] <- newClusNodes[1]
  newClusMRCAs <- c(newClusMRCAs, tail(newClusNodes, n = -1))
  
  newLogLik <- clusSplitMergeLogLik(AugTreePointer = currentValue$extPointer, clusMRCAsToSplitOrMerge = currentValue$paraValues$clusterNodeIndices[[clusNumber]], numOpenMP = numLikThreads, withinTransProbs = withinTransMatList, betweenTransProbs = betweenTransMatList, limProbs = limProbs)
  
  shortRecursive <- function(cMRCAs, cNumbers, clusInd, index = 1) {
    seqLabelsToChange <- currentPhylo$tip.label[phangorn::Descendants(currentPhylo, cMRCAs[index])[[1]]]
    clusInd[seqLabelsToChange] <- cNumbers[index]
    if (index <= length(cMRCAs)) {
      clusInd <- shortRecursive(cMRCAs, cNumbers, clusInd, index = index+1)
    }
    clusInd
  }
  newClusInd <- shortRecursive(newClusNodes, c(clusNumber, (length(currentValue$paraValues$clusterNodeIndices)+1):length(newClusMRCAs)), currentValue$paraValues$clusInd) ## I wonder if this is more efficient than apply + <<-...
  newCounts <- table(newClusInd)
  newLogPostProb <- newLogLik + clusIndLogPrior(clusInd = newClusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE) + dpois(max(newClusInd), lambda = poisRateNumClus, log = TRUE) ## Added the Poisson log-prob. to reflect a Poisson prior on the total number of clusters.
  
  tipAncestorsNew <- phangorn::Ancestors(currentPhylo, node = newClusMRCAs, type = "parent")
  ancestorsGroupsNew <- split(tipAncestorsNew, f = tipAncestorsNew)
  ancestorsGroupsPairsNew <- ancestorsGroupsNew[vapply(ancestorsGroupsNew, FUN = function(x) {length(x) > 1}, FUN.VALUE = c(Return = TRUE))]
  numPairsNew <- length(ancestorsGroupsPairsNew)
  numSplitsNew <- sum(newCounts > 1)
  transKernRatio <- numMoves/(numPairsNew + numSplitsNew)
  list(logLik = newLogLik, clusInd = newClusInd, counts = newCounts, logPostProb = newLogPostProb, transKernRatio = transKernRatio, clusterNodeIndices = newClusMRCAs)
}

.performMerge <- function(currentValue, clusMRCAsToMerge, DNAdataBin, limProbs, withinTransMatList, betweenTransMatList, numLikThreads, tipAncestors, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin, numMoves) {
    
  currentClusMRCAs <- currentValue$paraValues$clusterNodeIndices
  currentPhylo <- currentValue$paraValues$phylogeny
  
  newMRCA <- phangorn::mrca.phylo(currentPhylo, clusMRCAsToMerge)
  clusToMergeNumbers <- match(clusMRCAsToMerge, currentClusMRCAs)
  newClusMRCAs <- unique(replace(currentClusMRCAs, clusToMergeNumbers, newMRCA)) # No two clusters can have the same MRCA.
  
  newLogLik <- clusSplitMergeLogLik(AugTreePointer = currentValue$extPointer, clusMRCAsToSplitOrMerge = clusMRCAsToMerge, numOpenMP = numLikThreads, withinTransProbs = withinTransMatList, betweenTransProbs = betweenTransMatList, limProbs = limProbs)
  
  newClusInd <- replace(currentValue$paraValues$clusInd, which(currentValue$paraValues$clusInd %in% clusToMergeNumbers), min(clusToMergeNumbers)) ## New cluster takes the lowest index of the merged clusters, creating a gap.
  newCounts <- table(newClusInd)
  newLogPostProb <- newLogLik + clusIndLogPrior(clusInd = newClusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE) + dpois(max(newClusInd), lambda = poisRateNumClus, log = TRUE) ## Added the Poisson log-prob. to reflect a Poisson prior on the total number of clusters.
  tipAncestorsNew <- phangorn::Ancestors(currentPhylo, node = newClusMRCAs, type = "parent")
  ancestorsGroupsNew <- split(tipAncestorsNew, f = tipAncestorsNew)
  ancestorsGroupsPairsNew <- ancestorsGroupsNew[vapply(ancestorsGroupsNew, FUN = function(x) {length(x) > 1}, FUN.VALUE = c(Return = TRUE))]
  numPairsNew <- length(ancestorsGroupsPairsNew)
  numSplitsNew <- sum(newCounts > 1)
  transKernRatio <- numMoves/(numSplitsNew + numPairsNew)
  list(logLik = newLogLik, clusInd = newClusInd, counts = newCounts, logPostProb = newLogPostProb, transKernRatio = transKernRatio, clusterNodeIndices = newClusMRCAs)
}

.splitJoinClusterMove <- function(currentValue, DNAdataBin, limProbs, withinTransMatList, betweenTransMatList, numLikThreads, poisRateNumClus, shapeForAlpha, scaleForAlpha, alphaMin) {

    if (length(currentValue$paraValues$clusterNodeIndices) > 1) 
    {
      tipAncestors <- phangorn::Ancestors(currentValue$paraValues$phylogeny, node = currentValue$paraValues$clusterNodeIndices, type = "parent")
      ancestorsGroups <- split(currentValue$paraValues$clusterNodeIndices, f = tipAncestors)
      ancestorsGroupsPairs <- ancestorsGroups[vapply(ancestorsGroups, FUN = function(x) {length(x) > 1}, FUN.VALUE = logical(1))]
      numPairs <- length(ancestorsGroupsPairs)
    } 
    else 
    {
      numPairs <- 0
    }
    
    clustersToSplit <- names(currentValue$clusterCounts)[currentValue$clusterCounts > 1] 
    numSplits <- length(clustersToSplit)
   
    numMoves <- numPairs+numSplits
    if (numMoves == 0)
    {
      stop("No possible move. This is an error. Exiting... \n")
    }
    selectedMove <- sample(1:numMoves, size = 1)
    splitMove <- FALSE
    if (selectedMove > numPairs)
    {
      splitMove <- TRUE
    }
    
    if (splitMove) 
    {
      clusNumber <- as.numeric(clustersToSplit)[selectedMove - numPairs]
      splitMergeResult <- .performSplit(currentValue = currentValue, numMoves = numMoves, clusNumber = clusNumber, DNAdataBin = DNAdataBin, limProbs = limProbs, withinTransMatList = withinTransMatList, betweenTransMatList = betweenTransMatList, numLikThreads = numLikThreads, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
    } 
    else
    {
      pairToSelect <- selectedMove
      splitMergeResult <- .performMerge(currentValue = currentValue, clusMRCAsToMerge = ancestorsGroupsPairs[[pairToSelect]], DNAdataBin = DNAdataBin, limProbs = limProbs, withinTransMatList = withinTransMatList, betweenTransMatList = betweenTransMatList, numLikThreads = numLikThreads, tipAncestors = tipAncestors, poisRateNumClus = poisRateNumClus, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin, numMoves = numMoves)
    }
    MHratio <- splitMergeResult$transKernRatio*exp(splitMergeResult$logPostProb - currentValue$logPostProb)
    
    testUnif <- runif(1)
    if (testUnif < MHratio)
    {
      currentValue$logLik <- splitMergeResult$logLik
      currentValue$logPostProb <- splitMergeResult$logPostProb
      currentValue$paraValues$clusInd <- splitMergeResult$clusInd
      currentValue$clusterCounts <- splitMergeResult$counts
      currentValue$paraValues$clusterNodeIndices <- splitMergeResult$clusterNodeIndices
      
      if (!splitMove)
      {
        names(currentValue$clusterCounts) <- as.character(.relabel(as.numeric(names(currentValue$clusterCounts))))
        currentValue$paraValues$clusInd <- .relabel(currentValue$paraValues$clusInd)
      }
    } else {
      RestorePreviousConfig(currentValue$extPointer, edgeMat = matrix(0,1,1), NNImove = FALSE, withinMatListIndex = currentValue$paraValues$withinMatListIndex, betweenMatListIndex = currentValue$paraValues$betweenMatListIndex, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, splitMergeMove = TRUE) # Function tells the C++ AugTree object to restore the old configuration, also stored in the object.
    }
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
  
  if (betweenBool)
  {
    newLogLik <- newBetweenTransProbsLogLik(AugTreePointer = currentValue$extPointer, newBetweenTransProbs = betweenTransMatList, numOpenMP = numLikThreads, newBetweenMatListIndex = newIndex, withinTransProbs = withinTransMatList, limProbs = limProbs)
  }
  else
  {
    newLogLik <- newWithinTransProbsLogLik(AugTreePointer = currentValue$extPointer, newWithinTransProbs = withinTransMatList, numOpenMP = numLikThreads, newWithinMatListIndex = newIndex, betweenTransProbs = betweenTransMatList, limProbs = limProbs)
  }
 
  MHratio <- exp(newLogLik - currentValue$logLik)
  if (runif(1) < MHratio) {
    if (betweenBool) {
      currentValue$paraValues$betweenMatListIndex <- newIndex
    } else{
      currentValue$paraValues$withinMatListIndex <- newIndex
    }
    currentValue$logPostProb <- newLogLik +  (currentValue$logPostProb - currentValue$logLik)
    currentValue$logLik <- newLogLik
  } else {
    RestorePreviousConfig(currentValue$extPointer, edgeMat = matrix(0,1,1), NNImove = FALSE, withinMatListIndex = currentValue$paraValues$withinMatListIndex, betweenMatListIndex = currentValue$paraValues$betweenMatListIndex, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, splitMergeMove = FALSE) # Function tells the C++ AugTree object to restore the old configuration, also stored in the object.
  }
  currentValue
}

.DMphyClusCore <- function(nIter, startingValues, limProbs, shapeForAlpha, scaleForAlpha, numMovesNNIbetween, numMovesNNIwithin, numLikThreads, alignment, poisRateNumClus, clusPhyloUpdateProp, numSplitMergeMoves, alphaMin, withinTransMatAll, betweenTransMatAll, saveFrequency, intermediateDirectory = NULL, initialParaValues = NULL, maxMapSize = NULL, cullProportion = NULL, togglePhyloUpdates = TRUE) {

  if (is.matrix(withinTransMatAll[[1]])) {
    numGammaCat <- length(withinTransMatAll) ## We only have one set of substitution rate matrices.
  } else{
    numGammaCat <- length(withinTransMatAll[[1]])
  }
  DNAdataBin <- getConvertedAlignment(alignmentAlphaMat = alignment, equivVector = names(limProbs))
  names(DNAdataBin) <- rownames(alignment)
  
  if (is.null(initialParaValues))
  {
    currentValue <- .initCurrentValue(startingValues = startingValues, DNAdataBin = DNAdataBin, withinTransMatAll = withinTransMatAll, betweenTransMatAll = betweenTransMatAll, limProbs = limProbs, numLikThreads = numLikThreads, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin, poisRateNumClus = poisRateNumClus)
  } 
  else
  {
    currentValue <- initializeFromParameters(initialParaValues, withinTransMatAll, betweenTransMatAll, limProbs, numLikThreads, DNAdataBin = DNAdataBin, shapeForAlpha = shapeForAlpha, scaleForAlpha = scaleForAlpha, alphaMin = alphaMin)
  }

  cat("Launching chain... \n \n")

  progress <- txtProgressBar(style = 3)
  setTxtProgressBar(pb = progress, value = 0)
  onePerc <- nIter/100
  longOut <- lapply(1:nIter, FUN = function(x) {
    if ((x %% onePerc) == 0) {
      setTxtProgressBar(pb = progress, value = x/nIter)
    } else{}

    currentValue <<- .performStepPhylo(currentValue = currentValue, limProbs = limProbs, shapePriorAlpha = shapeForAlpha, scalePriorAlpha = scaleForAlpha, withinTransMatAll = withinTransMatAll, betweenTransMatAll = betweenTransMatAll, currentIter = x, numMovesNNIbetween = numMovesNNIbetween, numMovesNNIwithin = numMovesNNIwithin, numLikThreads = numLikThreads, DNAdataBin = DNAdataBin, poisRateNumClus = poisRateNumClus, clusPhyloUpdateProp = clusPhyloUpdateProp, alphaMin = alphaMin, numSplitMergeMoves = numSplitMergeMoves, maxMapSize = maxMapSize, cullProportion = cullProportion, togglePhyloUpdates = togglePhyloUpdates)

    paraVec <- currentValue$paraValues
    output <- c(paraVec, list(logPostProb = currentValue$logPostProb), list(logLik = currentValue$logLik))
    if (!is.null(intermediateDirectory))
    {
      if ((x %% saveFrequency) == 0) {
        filename <- paste(intermediateDirectory, "/simResults", digest::digest(list(.performStepPhylo, x)), "iter", x, ".Rdata", sep = "")
        save(output, file = filename, compress = TRUE)
      }
    }
    output
  })
  setTxtProgressBar(pb = progress, value = 1)
  close(con = progress)
  cat("\n Chain complete. \n\n\n")
  finalDeallocate(currentValue$extPointer)
  rm(currentValue) # Probably not needed... Once the function is done running, currentValue should go out of scope, and the destructor for AugTree should get called.
  longOut
}

initializeFromParameters <- function(initialParaValues, withinTransMatAll, betweenTransMatAll, limProbs, numLikThreads, DNAdataBin, alphaMin, shapeForAlpha, scaleForAlpha, poisRateNumClus)
{
  currentValue$paraValues <- initialParaValues
  currentValue$logLik <- logLikCpp(edgeMat = currentValue$paraValues$phylogeny$edge, limProbsVec = limProbs, withinTransMatList = withinTransMatAll[[currentValue$paraValues$withinMatListIndex]], betweenTransMatList = betweenTransMatAll[[currentValue$paraValues$betweenMatListIndex]], numOpenMP = numLikThreads, alignmentBin = DNAdataBin, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, numTips = ape::Ntip(currentValue$paraValues$phylogeny), numLoci = length(DNAdataBin[[1]]), withinMatListIndex = currentValue$paraValues$withinMatListIndex, betweenMatListIndex = currentValue$paraValues$betweenMatListIndex)
  currentValue$logPostProb <- currentValue$logLik + clusIndLogPrior(clusInd = currentValue$paraValues$clusInd, alpha = currentValue$paraValues$alpha) + dgamma(currentValue$paraValues$alpha - alphaMin, shape = shapeForAlpha, scale = scaleForAlpha, log = TRUE) + dpois(max(currentValue$paraValues$clusInd), lambda = poisRateNumClus, log = TRUE)
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

.updateBetweenPhylo <- function(currentValue, limProbs, withinTransMatList, betweenTransMatList, numMovesNNI, numLikThreads, DNAdataBin) {
    
  updatedPhyloAndLogLik <- betweenClusNNIlogLik(AugTreePointer = currentValue$extPointer, numMovesNNI = numMovesNNI, numOpenMP = numLikThreads, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, withinTransProbs = withinTransMatList, betweenTransProbs = betweenTransMatList, limProbs = limProbs)
  updatedLogLik <- updatedPhyloAndLogLik$logLik
  newPhylo <- list(edge = updatedPhyloAndLogLik$edge, tip.label = currentValue$paraValues$phylogeny$tip.label, edge.length = NULL, Nnode = ape::Nnode(currentValue$paraValues$phylogeny))
  class(newPhylo) <- "phylo"
  
  MHratio <- exp(updatedLogLik - currentValue$logLik) ## Prior doesn't change...
  if (runif(1) < MHratio) {
    currentValue$logPostProb <- currentValue$logPostProb - currentValue$logLik + updatedLogLik
    currentValue$logLik <- updatedLogLik
    currentValue$paraValues$phylogeny <- newPhylo
  } else {
    RestorePreviousConfig(currentValue$extPointer, edgeMat = currentValue$paraValues$phylogeny$edge, NNImove = TRUE, withinMatListIndex = currentValue$paraValues$withinMatListIndex, betweenMatListIndex = currentValue$paraValues$betweenMatListIndex, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, splitMergeMove = FALSE) # Function tells the C++ AugTree object to restore the old configuration, also stored in the object.
  }
  currentValue
}

getNNIbetweenPhylo <- function(phylogeny, clusterMRCAs, numMovesNNI) {
  phylogeny$node.label <- rep("0", Nnode(phylogeny))
  clusterMRCAsNodes <- clusterMRCAs[clusterMRCAs>ape::Ntip(phylogeny)]
  phylogeny$node.label[clusterMRCAsNodes - ape::Ntip(phylogeny)] <- as.character(clusterMRCAsNodes)
  internalPhylo <- phylogeny
  clusterPhylos <- lapply(clusterMRCAsNodes, FUN = function(x) {
    internalPhylo <<- drop.tip(internalPhylo, tip = phylogeny$tip.label[phangorn::Descendants(phylogeny, x)[[1]]], trim.internal = FALSE)
    ape::extract.clade(phylogeny, node = x)
  })
  names(clusterPhylos) <- as.character(clusterMRCAsNodes)
  neighbourTree <- phangorn::rNNI(internalPhylo, moves = numMovesNNI) # Make this a C++ function eventually...
  newPhylo <- internalPhylo
  lapply(names(clusterPhylos), FUN = function(x) {
    newPhylo <<- ape::bind.tree(x = newPhylo, y = clusterPhylos[[x]], where = match(x, newPhylo$tip.label))
  })
  newPhylo
}

.updateClusterPhylos <- function(currentValue, limProbs, betweenTransMatList, withinTransMatList, numMovesNNI = 1, numLikThreads, DNAdataBin, clusPhyloUpdateProp = 1) 
{
    if (clusPhyloUpdateProp < 1) 
    {
      treesSelectNum <- ceiling(length(currentValue$paraValues$clusterNodeIndices)*clusPhyloUpdateProp)
      selectNodeIndices <- sample(currentValue$paraValues$clusterNodeIndices, size = treesSelectNum, replace = FALSE) ## Uniformly selecting cluster-specific phylogenies for updates does not affect the transition kernel ratio. Also, I doubt it will greatly affect autocorrelation in the chain for cluster assignment indices.
    }
    else
    {
      selectNodeIndices <- currentValue$paraValues$clusterNodeIndices
    }
    originalPhylo <- currentValue$paraValues$phylogeny # Meant to maintain the link between the numeric cluster MRCAs and corresponding tips.
    performMHupdate <- function(clusterMRCA) 
    {
      nodeDescendants <- originalPhylo$tip.label[phangorn::Descendants(x = originalPhylo, node = clusterMRCA)[[1]]]
      if (length(nodeDescendants) < 3) 
      {
        return(NULL) # NNI moves cannot be performed for phylogenies with 1 or two tips.
      }
      
      newLogLikAndPhylo <- withinClusNNIlogLik(AugTreePointer = currentValue$extPointer, MRCAofClusForNNI = clusterMRCA, numMovesNNI = numMovesNNI, numOpenMP = numLikThreads, withinTransProbs = withinTransMatList, betweenTransProbs = betweenTransMatList, limProbs = limProbs)
      newLogLik <- newLogLikAndPhylo$logLik
      newEdge <- newLogLikAndPhylo$edge
      newBigPhylo <- list(edge = newEdge, tip.label = currentValue$paraValues$phylogeny$tip.label, edge.length = NULL, Nnode = ape::Nnode(currentValue$paraValues$phylogeny))
      class(newBigPhylo) <- "phylo"
      newClusterMRCAs <- currentValue$paraValues$clusterNodeIndices # A NNI move should not change the numbering of the cluster root nodes. 
      
      MHratio <- exp(newLogLik - currentValue$logLik)
      
      if (runif(1) < MHratio) 
      {
        currentValue$paraValues$clusterNodeIndices <<- newClusterMRCAs
        currentValue$logPostProb <<- currentValue$logPostProb - currentValue$logLik + newLogLik
        currentValue$logLik <<- newLogLik
        currentValue$paraValues$phylogeny <<- newBigPhylo
      } 
      else
      {
        RestorePreviousConfig(currentValue$extPointer, edgeMat = currentValue$paraValues$phylogeny$edge, NNImove = TRUE, withinMatListIndex = currentValue$paraValues$withinMatListIndex, betweenMatListIndex = currentValue$paraValues$betweenMatListIndex, clusterMRCAs = currentValue$paraValues$clusterNodeIndices, splitMergeMove = FALSE)
      }
      NULL
    }
    lapply(selectNodeIndices, FUN = performMHupdate)
    currentValue
}

.checkInput <- function(startingValues, Qmatrix, alignment, limProbs, shiftForAlpha) {
    ## Check if all necessary starting values are defined...
    test <- all(c("alpha", "clusInd", "phylogeny") %in% names(startingValues))
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
    if (!is.null(Qmatrix)) {
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