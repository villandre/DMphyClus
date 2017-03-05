#include "Forest.h"
#include "InputNode.h"
#include "IntermediateNode.h"
#include "TreeNode.h"
#include <unordered_map>
#include "helper.h"

using namespace Rcpp ;
using namespace arma ;


Forest::Forest(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, std::vector<std::vector<uvec>> * alignmentBinPoint, const List & withinTransProbMatList, const List & betweenTransProbMatList, const NumericVector & limProbs, const uint numTips, const uint numLoci, solutionDictionaryType & solutionDictionary, const uint withinMatListIndex, const uint betweenMatListIndex)
{
  _numLoci = numLoci ;
  _numTips = numTips ;
  _numRateCats = withinTransProbMatList.size() ;
  _solutionDictionary = solutionDictionary ;
  _forest.reserve(alignmentBinPoint->size()*withinTransProbMatList.size()) ;
  _withinTransProbMatVec = as<std::vector<mat>>(withinTransProbMatList) ;
  _betweenTransProbMatVec = as<std::vector<mat>>(betweenTransProbMatList) ;
  umat edgeMatrixRecast = as<umat>(edgeMatrix) ;
  uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
  _limProbs = as<vec>(limProbs) ;
  _withinMatListIndex = withinMatListIndex ;
  _betweenMatListIndex = betweenMatListIndex ;
  
  _randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ; // This is the random number generator. It's initialized when the Forest is built, and the seed is 0 by default.
  _alignmentBinReference = alignmentBinPoint ;
 
  for (auto & i : *alignmentBinPoint) // Iterating on loci...
  {
    uint locusRateIndex = 0 ;
    
    for (uint rateCategIndex = 0 ; rateCategIndex < _numRateCats ; rateCategIndex++)
    {
      AugTree * LocusRateAugTree = new AugTree(edgeMatrixRecast, clusterMRCAsRecast, &i, rateCategIndex, _solutionDictionary, _numTips) ;
      LocusRateAugTree->ComputeKeys(LocusRateAugTree->GetVertexVector().at(_numTips), solutionDictionary, _withinMatListIndex, _betweenMatListIndex) ; // We start obtaining keys at the root.
      _forest.push_back(LocusRateAugTree) ;
    }
  }
}

Forest::Forest(const IntegerMatrix & edgeMatrix, const vec & limProbs, uint numRateCats, uint numLoci, uint numTips, gsl_rng * ranNumGenerator, solutionDictionaryType solutionDictionary, std::vector<std::vector<uvec>> * alignmentBinPoint, const uint withinMatListIndex, const uint betweenMatListIndex)
{ 
  _numRateCats = numRateCats ;
  _numLoci = numLoci ;
  _numTips = numTips ;
  _randomNumGenerator = ranNumGenerator ;
  _forest.reserve(numLoci*numRateCats) ;
  umat edgeMatrixRecast = as<umat>(edgeMatrix) ;
  _solutionDictionary = solutionDictionary ;
  _alignmentBinReference = alignmentBinPoint ;
  _withinMatListIndex = withinMatListIndex ;
  _betweenMatListIndex = betweenMatListIndex ;
  
  for (uint i = 0; i < numLoci; i++) // Iterating on loci...
  {
    for (uint j = 0 ; j < numRateCats ; j++)
    {
      AugTree * LocusRateAugTree = new AugTree(edgeMatrixRecast, j, _numTips) ;
      _forest.push_back(LocusRateAugTree) ; 
    }
  }
}

void Forest::ComputeLoglik()
{
  uint rateCategIndex = 0 ;
  //#pragma omp parallel for 
  for (std::vector<AugTree *>::iterator forestIter = _forest.begin(); forestIter < _forest.end(); forestIter++) 
  {
    (*forestIter)->SolveRoot(_solutionDictionary, _withinTransProbMatVec.at(rateCategIndex), _betweenTransProbMatVec.at(rateCategIndex), _limProbs, _numTips) ;
    rateCategIndex = littleCycle(rateCategIndex+1, _withinTransProbMatVec.size()) ;
  }
  // Now, we must average likelihoods across rate categories for each locus, log the output, and sum the resulting logs.
  vec rateAveragedLogLiks(_numLoci) ;
  vec likAcrossRatesLoci(_forest.size()) ;
  vec exponentVec(_forest.size()) ;
  
  std::transform(_forest.begin(), _forest.end(), likAcrossRatesLoci.begin(), [] (AugTree * myTree) {return myTree->GetLikelihood() ;}) ;
  std::transform(_forest.begin(), _forest.end(), exponentVec.begin(), [] (AugTree * myTree) {return myTree->GetExponentContainer() ;}) ;
  
  for (uint i = 0; i < rateAveragedLogLiks.size(); i++)
  {
    double maxExponent = max(exponentVec.rows(_numRateCats*i, _numRateCats*(i+1) - 1)) ;
    exponentVec.rows(_numRateCats*i, _numRateCats*(i+1) - 1) -= maxExponent ;
    rateAveragedLogLiks[i] = log(mean(likAcrossRatesLoci.rows(_numRateCats*i, _numRateCats*(i+1) - 1)%exp(exponentVec.rows(_numRateCats*i, _numRateCats*(i+1) - 1)))) + maxExponent;
  }
  
  _loglik = sum(rateAveragedLogLiks) ;
}

void Forest::HandleSplit(uint clusMRCAtoSplit)
{
  for (auto & augtree : _forest)
  {
    augtree->GetVertexVector().at(clusMRCAtoSplit - 1)->InvalidateSolution() ;
    
    for (auto & childNode : augtree->GetVertexVector().at(clusMRCAtoSplit - 1)->GetChildren())
    {
      childNode->SetWithinParentBranch(false) ;
    }
    augtree->ComputeKeys(augtree->GetVertexVector().at(_numTips), _solutionDictionary, _withinMatListIndex, _betweenMatListIndex) ;
    augtree->PatternLookup(_solutionDictionary, augtree->GetVertexVector().at(_numTips)) ;
  }
}

void Forest::HandleMerge(uvec & clusMRCAstoMerge)
{
  for (auto & augtree : _forest)
  {
    augtree->GetVertexVector().at(clusMRCAstoMerge.at(0) - 1)->GetParent()->InvalidateSolution() ; // Elements of clusMRCAsToMerge should all have the same parent to allow a merge to occur.
    
    for (auto & oldClusterMRCA : clusMRCAstoMerge)
    {
      augtree->GetVertexVector().at(oldClusterMRCA - 1)->SetWithinParentBranch(true) ;
    }
    augtree->ComputeKeys(augtree->GetVertexVector().at(_numTips), _solutionDictionary, _withinMatListIndex, _betweenMatListIndex) ;
    augtree->PatternLookup(_solutionDictionary, augtree->GetVertexVector().at(_numTips)) ;
  }
}

void Forest::InputForestElements(Forest * originForest)
{
  _withinTransProbMatVec = originForest->GetWithinTransProbMatVec() ;
  _betweenTransProbMatVec = originForest->GetBetweenTransProbMatVec() ;
  _withinMatListIndex = originForest->GetWithinMatListIndex() ;
  _betweenMatListIndex = originForest->GetBetweenMatListIndex() ;
  uint originAugTreeIndex = 0 ;
  for (auto & i : _forest)  
  {
    i->CopyAugTreeNonPointer(originForest->GetForest().at(originAugTreeIndex)) ;
    originAugTreeIndex++ ;
  }
}

void Forest::InvalidateBetweenSolutions()
{
  for (auto & augtree : _forest)
  {
    augtree->CheckAndInvalidateBetweenRecursive(augtree->GetVertexVector().at(_numTips)) ;
    augtree->ComputeKeys(augtree->GetVertexVector().at(_numTips), _solutionDictionary, _withinMatListIndex, _betweenMatListIndex) ;
    augtree->PatternLookup(_solutionDictionary, augtree->GetVertexVector().at(_numTips)) ;
  }
}

void Forest::InvalidateAllSolutions()
{
  for (auto & augtree : _forest)
  {
    augtree->InvalidateAll() ;
    augtree->ComputeKeys(augtree->GetVertexVector().at(_numTips), _solutionDictionary, _withinMatListIndex, _betweenMatListIndex) ;
    augtree->PatternLookup(_solutionDictionary, augtree->GetVertexVector().at(_numTips)) ;
  }
}

void Forest::RearrangeNNI(const uint vertexId1, const uint vertexId2)
{
  for (auto & i : _forest)
  {
    i->RearrangeTreeNNI(vertexId1, vertexId2, _solutionDictionary) ;
    i->ComputeKeys(i->GetVertexVector().at(_numTips), _solutionDictionary, _withinMatListIndex, _betweenMatListIndex) ;
    i->PatternLookup(_solutionDictionary, i->GetVertexVector().at(_numTips)) ;
  }
}

void Forest::RebuildTrees(const umat & edgeMat)
{
  for (auto & i : _forest)
  {
    i->BuildTreeNoAssign(edgeMat) ;
  }
}