% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMphyClusInterface.R
\name{DMphyClusChain}
\alias{DMphyClusChain}
\title{Cluster an alignment with DM-PhyClus}
\usage{
DMphyClusChain(numIters, numLikThreads = 1, numMovesNNIbetween = 1,
  numMovesNNIwithin = 1, alignment, startingValues,
  numSamplesForTransMat = 1e+05, coefVarForTransMat = 1,
  meanBetweenBranchVec, meanWithinBranchVec, limProbs,
  clusPhyloUpdateProp = 1, numSplitMergeMoves = 1, numGammaCat = 3,
  discGammaPar = NULL, Qmatrix = NULL, shapeForAlpha, scaleForAlpha,
  shiftForAlpha = 0, poisRateNumClus, betweenClusTransMatList = NULL,
  withinClusTransMatList = NULL, intermediateDirectory = NULL,
  saveFrequency = 20, initialParaValues = NULL, maxMemory = NULL,
  cullProportion = NULL, numClusters = NULL)
}
\arguments{
\item{numIters}{number of iterations for the MCMC sampler}

\item{numLikThreads}{number of openMP threads for log-likelihood evaluations}

\item{numMovesNNIbetween}{number of nearest-neighbour interchange moves
used in proposing transitions in the space of between-cluster phylogenies}

\item{numMovesNNIwithin}{number of nearest-neighbour interchange moves used
in proposing transitions in the space of within-cluster phylogenies}

\item{alignment}{matrix of characters. When using DNA data, all nucleotides
and ambiguities should be coded with small letters, following the IUPAC standard}

\item{startingValues}{list of starting values for the model parameters: elements
named phylogeny, clusInd, and alpha are required}

\item{numSamplesForTransMat}{number of random values used in Monte Carlo estimation
of the transition probability matrices}

\item{meanBetweenBranchVec}{vector of positive numerical values giving a range of
potential mean branch lengths in the between-cluster phylogeny}

\item{meanWithinBranchVec}{vector of positive numerical values giving a range of
potential mean branch lengths in the within-cluster phylogenies}

\item{limProbs}{vector of numerical values, with named elements, giving the limiting
probabilities for all states}

\item{clusPhyloUpdateProp}{number between 0 and 1 indicating the proportion of
within-cluster phylogenies to update in each MCMC iteration}

\item{numSplitMergeMoves}{the number of times the algorithm should try to split or
merge clusters in each MCMC iteration}

\item{numGammaCat}{number of among-loci substitution rate variation categories}

\item{discGammaPar}{value of the the discrete gamma distribution parameter tuning
among-sites rate variation}

\item{Qmatrix}{substitution rate matrix with rows and columns names corresponding
to states in the alignment}

\item{shapeForAlpha}{shape parameter of the gamma distribution used as a prior for
the Dirichlet concentration parameter}

\item{scaleForAlpha}{scale parameter of the gamma distribution used as a prior for
the Dirichlet concentration parameter}

\item{shiftForAlpha}{minimum value the Dirichlet concentration parameter can take}

\item{poisRateNumClus}{Poisson parameter for the number of clusters distribution;
if left unspecified, it is set equal to the number of clusters in the startingValues}

\item{betweenClusTransMatList}{list of lists of potential transition rate matrices
along branches in the supporting phylogeny. See details}

\item{withinClusTransMatList}{Same as betweenClusTransMatList, but for the
within-cluster phylogenies}

\item{intermediateDirectory}{Directory where intermediate results will be saved. If left unspecified, intermediate results will not be saved.}

\item{saveFrequency}{Defaults to 20. Determines the frequency at which intermediate outputs will be saved. Does not apply if intermediateDirectory is left unspecified.}

\item{maxMemory}{maximum memory size (in megs).}

\item{cullProportion}{the proportion of entries that must be removed when the solution dictionary is culled.}

\item{numClusters}{(experimental) the number of clusters}
}
\value{
A list with two components:
\itemize{
 \item{chain:} {list with each element itself a list giving the sampled parameter values.
betweenTransMatListIndex and withinTransMatListIndex correspond to the index of the assumed mean
branch length in meanBetweenBranchVec and meanWithinBranchVec, respectively}
 \item{MAPestimate:}{ vector giving the maximum posterior probability cluster membership indices estimate}
}
}
\description{
Implementation of the DM-PhyClus method described in Villandr\'{e} et al. 2016.
}
\details{
It is possible to supply transition probability matrices directly to the function
by specifying values for betweenClusTransMatList and withinClusTransMatList. Since estimating
those matrices might take a long time, it is preferable to perform the computations only once
using the \code{\link{outputTransMatList}} function, save the result to a file and re-use the
matrices whenever required. When those matrices are provided, no value need be specified for
Qmatrix, numGammaCat, discGammaPar, meanBetweenBranchVec, and meanWithinBranchVec.
}
\examples{
\dontrun{
INPUT_AN_EXAMPLE()
}
}
