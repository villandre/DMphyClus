#include "AugTree.h"
#include "InputNode.h"
#include "IntermediateNode.h"
#include "TreeNode.h"
#include <unordered_map>
#include "helper.h"

using namespace Rcpp ;
using namespace arma ;

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}


AugTree::AugTree(const umat & edgeMatrix, const uvec & clusterMRCAs, std::vector<uvec> & alignmentBinOneLocusOneRate, const Col<double> & limProbs, const uint numTips, const uint rateCategIndex, solutionDictionaryType & solutionDictionary)
{ 
  _exponentContainer = 0 ;
  _numTips = numTips ;
  _limProbs = limProbs ;
  _rateCateg = rateCategIndex ;
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
  InitializeVertices(alignmentBinOneLocusOneRate) ;
  AssociateTransProbMatrices(clusterMRCAs) ; 
//  ComputeKeys(_vertexVector[_numTips], solutionDictionary) ; // We start obtaining keys at the root.
}

AugTree::AugTree(const umat & edgeMatrix, const vec & limProbs, const uint numTips, const uint rateCategIndex)
{
  _exponentContainer = 0 ;
  _numTips = numTips ;
  _rateCateg = rateCategIndex ;
  _limProbs = limProbs ;
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
}

void AugTree::CopyAugTreeNonPointer(AugTree * sourceAugTree) 
{
  _exponentContainer = sourceAugTree->GetExponentContainer() ;
  uint sourceVertexIndex = 0 ; // Defining an iterator on sourceAugTree->GetVertexVector() cannot be done without copying it, which I can't explain for now.
  for (auto & i : _vertexVector)
  {
    i->EnterCommonInfo(sourceAugTree->GetVertexVector().at(sourceVertexIndex)) ;
    i->EnterSolution(sourceAugTree->GetVertexVector().at(sourceVertexIndex)) ;
    sourceVertexIndex++ ;
  }
}

void AugTree::AssociateTransProbMatrices(const uvec & clusterMRCAs) 
{
  // By default all nodes are considered between clusters.
  for (auto & i : _vertexVector)
  {
    i->SetWithinParentBranchAndRateCateg(false, _rateCateg) ;
  }
  
  for (auto & i : clusterMRCAs)
  {
    if (i > _numTips) { // Again, clusterMRCAs is based on the R convention, hence >, and not >=.
      for (auto & j : _vertexVector[i-1]->GetChildren()) {
        BindMatrix(j, true) ; 
      }  
    }
  }
}

void AugTree::BindMatrix(TreeNode * vertex, const bool withinCluster)
{
  vertex->SetWithinParentBranch(withinCluster) ;
  if (!(vertex->GetChildren()[0] == NULL)) { // A null pointer indicates that we've reached an input node.
    for (auto & i : vertex->GetChildren())
    {
      BindMatrix(i, withinCluster) ;
    }
  }
}

void AugTree::BuildTree(umat & edgeMatrix)
{
  _vertexVector.reserve(edgeMatrix.n_rows + 1) ;

  // We create the tips. Note that tip 1 should correspond to vertex 1 in the original (the one in the phylo object) edgeMatrix

  for (uint i = 0; i < _numTips; i++) {
    InputNode * newNode = new InputNode{};
    _vertexVector.push_back(newNode) ; 
  } ;

  // We add the internal nodes
  for (uint i = 0 ; i < edgeMatrix.n_rows - _numTips + 1; i++) {
    IntermediateNode * newNode = new IntermediateNode{};
    _vertexVector.push_back(newNode) ; 
  } ;
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _vertexVector.size(); i++) {
    _vertexVector[i]->SetId(i) ; 
  } ;
  // The vertices are all disjoint, the next loop defines their relationships
  // The iterator follows columns.
  for (umat::iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end()-edgeMatrix.n_rows; iter++)
  {
    _vertexVector[*iter]->AddChild(_vertexVector[*(iter+edgeMatrix.n_rows)]) ;
    _vertexVector[*(iter+edgeMatrix.n_rows)]->SetParent(_vertexVector[*iter]) ;
  }
}

void AugTree::ComputeKeys(TreeNode * vertex, solutionDictionaryType & solutionDictionary)
{
  if (!vertex->IsKeyDefined())
  {
    if (!vertex->CanFindKey()) // If a vertex can be solved, then its children have patterns assigned to them.
    {
      for (auto & i : vertex->GetChildren())
      {
        ComputeKeys(i, solutionDictionary) ;
      }
    }
    vertex->DeriveKey(solutionDictionary) ;
  }
}

void AugTree::SolveRoot(solutionDictionaryType & solutionDictionary, const mat & withinTransProbMat, const mat & betweenTransProbMat)
{
  //ComputeKeys(_vertexVector[_numTips], solutionDictionary) ;
  //PatternLookup(solutionDictionary, _vertexVector[_numTips]) ;
  TrySolve(_vertexVector[_numTips], solutionDictionary, withinTransProbMat, betweenTransProbMat) ;
  _likelihoodProp = dot(_vertexVector[_numTips]->GetSolution(), _limProbs) ;
}

void AugTree::InitializeVertices(std::vector<uvec> & alignmentBinOneLocusOneRate)
{ 
  std::vector<TreeNode *>::iterator treeIter = _vertexVector.begin() ;
  for (auto & i : alignmentBinOneLocusOneRate)
  {
    (*treeIter)->SetInput(&i) ; 
    //(*treeIter)->ToggleSolved() ;
    treeIter++ ;
  }
}

void AugTree::PatternLookup(solutionDictionaryType & solutionDictionary, TreeNode * currentNode) 
{
  if (!currentNode->IsSolved() && !solutionDictionary->empty())
  { // If the node is already solved, no need to update it with a stored pattern.
    if (solutionDictionary->count(currentNode->GetDictionaryKey()) == 0)
    {
      for (auto & i : currentNode->GetChildren())
      {
        PatternLookup(solutionDictionary, i) ;
      }
    }
    else
    {
      currentNode->SetSolution((*solutionDictionary)[currentNode->GetDictionaryKey()]) ;
      currentNode->SetSolved(true) ;
    }
  }
}

void AugTree::TrySolve(TreeNode * vertex, solutionDictionaryType & solutionDictionary, const mat & withinTransProbMat, const mat & betweenTransProbMat)
{
  if (!(vertex->IsSolved())) // Could be solved because of the pattern lookup.
  {
    if (!vertex->CanSolve())
    {
      for (auto & i : vertex->GetChildren())
      {
        TrySolve(i, solutionDictionary, withinTransProbMat, betweenTransProbMat) ;
      }
    }
    if (vertex->GetChildren().at(0)->GetWithinParentBranch()) 
    {  // This junction is within a cluster.
      vertex->ComputeSolution(solutionDictionary, withinTransProbMat, &_exponentContainer) ;
    }
    else
    {
      vertex->ComputeSolution(solutionDictionary, betweenTransProbMat, &_exponentContainer) ;
    }
  }
}

void AugTree::InvalidateAll() // Assumes that the tree starts fully solved.
{
  for (auto & i : _vertexVector)
  {
    i->SetSolved(false) ;
    i->MarkKeyUndefined() ;
  }
}

std::vector<uint> AugTree::GetNNIverticesWithin(TreeNode * originVertex)
{
  // We look at grandchildren. If the vertex has at least one grandchild, a NNI move from this vertex is possible.
  std::vector<uint> NNIvertices ;
  GetNNIverticesInternalWithin(originVertex, &NNIvertices) ;
  return NNIvertices ; 
}

void AugTree::GetNNIverticesInternalWithin(TreeNode * currentVertex, std::vector<uint> * NNIvertices) 
{
  bool movePossible = false ;
  if (currentVertex->GetChildren().at(0) != NULL) 
  {
    for (auto & i : currentVertex->GetChildren())
    {
      movePossible = movePossible || (i->GetChildren().at(0) != NULL) ;
      
      if (movePossible)
      { 
        break ; // As soon as this is true, logical "or" will not be able to yield false.
      }
    }
    if (movePossible)
    {
      NNIvertices->push_back(currentVertex->GetId()) ;
      for (auto & i : currentVertex->GetChildren())
      {
        GetNNIverticesInternalWithin(i, NNIvertices) ;
      }
    }
  }
}

std::vector<uint> AugTree::GetNNIverticesBetween(TreeNode * originVertex, uvec & clusterMRCAs)
{
  // We look at grandchildren. If the vertex has at least one grandchild, a NNI move from this vertex is possible.
  std::vector<uint> NNIvertices ;
  GetNNIverticesInternalBetween(originVertex, &NNIvertices, clusterMRCAs) ;
  return NNIvertices ; 
}

void AugTree::GetNNIverticesInternalBetween(TreeNode * currentVertex, std::vector<uint> * NNIvertices, uvec & clusterMRCAs) 
{
  bool movePossible = false ;
  
  if ((currentVertex->GetChildren().at(0) != NULL) && (std::find(clusterMRCAs.begin(), clusterMRCAs.end(), currentVertex->GetId()+1) == clusterMRCAs.end())) 
  {
    for (auto & i : currentVertex->GetChildren())
    {
      movePossible = (movePossible || ((i->GetChildren().at(0) != NULL) && (std::find(clusterMRCAs.begin(), clusterMRCAs.end(), i->GetId()+1) == clusterMRCAs.end()))) ;
      
      if (movePossible)
      { 
        break ; // As soon as this is true, logical "or" will not be able to yield false.
      }
    }
    if (movePossible)
    {
      NNIvertices->push_back(currentVertex->GetId()) ;
      for (auto & i : currentVertex->GetChildren())
      {
        GetNNIverticesInternalBetween(i, NNIvertices, clusterMRCAs) ;
      }
    }
  }
}


void AugTree::RearrangeTreeNNI(uint vertexId1, uint vertexId2) 
{
  _vertexVector.at(vertexId1)->GetParent()->RemoveChild(_vertexVector.at(vertexId1)) ;
  _vertexVector.at(vertexId1)->GetParent()->AddChild(_vertexVector.at(vertexId2)) ;
  
  _vertexVector.at(vertexId2)->GetParent()->RemoveChild(_vertexVector.at(vertexId2)) ;
  _vertexVector.at(vertexId2)->GetParent()->AddChild(_vertexVector.at(vertexId1)) ;
  
  TreeNode * oriParent1 = _vertexVector.at(vertexId1)->GetParent() ;
  TreeNode * oriParent2 = _vertexVector.at(vertexId2)->GetParent() ;
  _vertexVector.at(vertexId1)->SetParent(oriParent2) ;
  _vertexVector.at(vertexId2)->SetParent(oriParent1) ;
  
  _vertexVector.at(vertexId1)->GetParent()->InvalidateSolution() ;
  _vertexVector.at(vertexId2)->GetParent()->InvalidateSolution() ;
}

umat AugTree::BuildEdgeMatrix()
{
  umat edgeMatrix(_vertexVector.size()-1, 2, fill::zeros) ;
  uint lineNum = 0 ;
  AddEdgeRecursion(edgeMatrix, lineNum, _vertexVector.at(_numTips)) ;
  return edgeMatrix ;
}

void AugTree::AddEdgeRecursion(umat & matToUpdate, uint & lineNum, TreeNode * currentVertex)
{
  if (currentVertex->GetParent() != NULL)
  {
    matToUpdate.at(lineNum, 0) = currentVertex->GetParent()->GetId()+1 ; // Vertex numbering starts at 1 in R, instead of 0.
    matToUpdate.at(lineNum, 1) = currentVertex->GetId()+1 ;
    lineNum = lineNum + 1 ;
  }
  if (currentVertex->GetChildren().at(0) != NULL) 
  {
    for (auto & i : currentVertex->GetChildren())
    {
      AddEdgeRecursion(matToUpdate, lineNum, i) ;
    }
  }
}

Forest::Forest(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, std::vector<std::vector<uvec>> * alignmentBinPoint, const List & withinTransProbMatList, const List & betweenTransProbMatList, const NumericVector & limProbs, const uint numTips, const uint numLoci, solutionDictionaryType & solutionDictionary)
{
  _numLoci = numLoci ;
  _numRateCats = withinTransProbMatList.size() ;
  _solutionDictionary = solutionDictionary ;
  _forest.reserve(alignmentBinPoint->size()*withinTransProbMatList.size()) ;
  _withinTransProbMatVec = as<std::vector<mat>>(withinTransProbMatList) ;
  _betweenTransProbMatVec = as<std::vector<mat>>(betweenTransProbMatList) ;
  umat edgeMatrixRecast = as<umat>(edgeMatrix) ;
  uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
  Col<double> limProbsRecast = as<Col<double>>(limProbs) ;
  
  _randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ; // This is the random number generator. It's initialized when the Forest is built, and the seed is 0 by default.
  _alignmentBinReference = alignmentBinPoint ;
  
  for (auto & i : *alignmentBinPoint) // Iterating on loci...
  {
    uint locusRateIndex = 0 ;
    
    for (uint rateCategIndex = 0 ; rateCategIndex < _numRateCats ; rateCategIndex++)
    {
      AugTree * LocusRateAugTree = new AugTree(edgeMatrixRecast, clusterMRCAsRecast, i, limProbsRecast, numTips, rateCategIndex, _solutionDictionary) ;
      _forest.push_back(LocusRateAugTree) ;
    }
  }
}

Forest::Forest(const IntegerMatrix & edgeMatrix, const vec & limProbs, uint numRateCats, uint numLoci, uint numTips, gsl_rng * ranNumGenerator, solutionDictionaryType solutionDictionary, std::vector<std::vector<uvec>> * alignmentBinPoint)
{ 
  _numRateCats = numRateCats ;
  _numLoci = numLoci ;
  _randomNumGenerator = ranNumGenerator ;
  _forest.reserve(numLoci*numRateCats) ;
  umat edgeMatrixRecast = as<umat>(edgeMatrix) ;
  _solutionDictionary = solutionDictionary ;
  _alignmentBinReference = alignmentBinPoint ;
  
  for (uint i = 0; i < numLoci; i++) // Iterating on loci...
  {
    for (uint j = 0 ; j < numRateCats ; j++)
    {
      AugTree * LocusRateAugTree = new AugTree(edgeMatrixRecast, limProbs, numTips, j) ;
      _forest.push_back(LocusRateAugTree) ; 
    }
  }
}


void Forest::ComputeLoglik()
{
  uint rateCategIndex = 0 ;
  #pragma omp parallel for 
  for (std::vector<AugTree *>::iterator forestIter = _forest.begin(); forestIter < _forest.end(); forestIter++) // This syntax is compatible with openMP, unlike the more conventional 'for (auto & i : myVec')
  {
    (*forestIter)->SolveRoot(_solutionDictionary, _withinTransProbMatVec.at(rateCategIndex), _betweenTransProbMatVec.at(rateCategIndex)) ;
    rateCategIndex = littleCycle(rateCategIndex+1, _withinTransProbMatVec.size()) ;
  }
  
  // Now, we must average likelihoods across rate categories for each locus, log the output, and sum the resulting logs.
  Col<double> rateAveragedLogLiks(_numLoci) ;
  Col<double> likAcrossRatesLoci(_forest.size()) ;
  Col<double> exponentVec(_forest.size()) ;
  
  std::transform(_forest.begin(), _forest.end(), likAcrossRatesLoci.begin(), [] (AugTree * myTree) {return myTree->GetLikelihood() ;}) ;
  std::transform(_forest.begin(), _forest.end(), exponentVec.begin(), [] (AugTree * myTree) {return myTree->GetExponentContainer() ;}) ;
  
  for (uint i = 0; i < rateAveragedLogLiks.size(); i++)
  {
    double maxExponent = max(exponentVec.rows(_numRateCats*i, _numRateCats*(i+1) - 1)) ;
    exponentVec.rows(_numRateCats*i, _numRateCats*(i+1) - 1) -= maxExponent ;
    rateAveragedLogLiks[i] = log(mean(likAcrossRatesLoci.rows(_numRateCats*i, _numRateCats*(i+1) - 1)%exp(exponentVec.rows(_numRateCats*i, _numRateCats*(i+1) - 1)))) + maxExponent;
  }
  rateAveragedLogLiks.print("Log-liks per rate:") ;
  
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
  }
}

std::vector<uint> AugTree::GetTwoVerticesForNNI(gsl_rng * randomNumGenerator, TreeNode * subtreeRoot, uvec & clusterMRCAs)
{
  std::vector<uint> grandChildrenVec ;
  unsigned long int selectedChildIndex ;
  
  for (auto & i : subtreeRoot->GetChildren()) // This only works for bifurcating trees... For multifurcating trees, the two branches for NNI must also be selected.
  {
    if (i->GetChildren().at(0) == NULL || (std::find(clusterMRCAs.begin(), clusterMRCAs.end(), i->GetId() + 1) != clusterMRCAs.end()))
    {
      grandChildrenVec.push_back(i->GetId()) ;
    }
    else
    {
      selectedChildIndex = gsl_rng_uniform_int(randomNumGenerator, i->GetChildren().size()) ;
      grandChildrenVec.push_back(i->GetChildren().at(selectedChildIndex)->GetId()) ;
    }
  }
  return grandChildrenVec ;
}

void Forest::InputForestElements(XPtr<Forest> originForest)
{
  _withinTransProbMatVec = originForest->GetWithinTransProbMatVec() ;
  _betweenTransProbMatVec = originForest->GetBetweenTransProbMatVec() ;
  uint originAugTreeIndex = 0 ;
  for (auto & i : _forest)  {
    
    i->CopyAugTreeNonPointer(originForest->GetForest().at(originAugTreeIndex)) ;
    originAugTreeIndex++ ;
  }
}

void Forest::InvalidateBetweenSolutions()
{
  uint numTips = _forest.at(0)->GetNumTips() ;
  for (auto & augtree : _forest)
  {
    augtree->CheckAndInvalidateBetweenRecursive(augtree->GetVertexVector().at(numTips)) ;
  }
}

void AugTree::CheckAndInvalidateBetweenRecursive(TreeNode * currentVertex) 
{
  if (currentVertex->GetChildren().at(0) != NULL)
  {
    if (!(currentVertex->GetChildren().at(0)->GetWithinParentBranch()))
    {
      currentVertex->SetSolved(false) ;
      currentVertex->MarkKeyUndefined() ;
      
      for (auto & child : currentVertex->GetChildren())
      {
        CheckAndInvalidateBetweenRecursive(child) ;
      }
    }
  }
}

void Forest::InvalidateAllSolutions()
{
  for (auto & augtree : _forest)
  {
    augtree->InvalidateAll() ;
  }
}