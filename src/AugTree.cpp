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


AugTree::AugTree(const umat & edgeMatrix, const uvec & clusterMRCAs, const mat & withinTransProbMatrix, const mat & betweenTransProbMatrix, const std::vector<uvec> & alignmentBinOneLocusOneRate, const Col<double> & limProbs, const uint numTips, const uint rateCategIndex, solutionDictionaryType & solutionDictionary)
{
  _numTips = numTips ;
  _limProbs = limProbs ;
  _rateCateg = rateCategIndex ;
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
  InitializeVertices(alignmentBinOneLocusOneRate) ;
  AssociateTransProbMatrices(clusterMRCAs, withinTransProbMatrix, betweenTransProbMatrix) ;
  //ComputeKeys(_vertexVector[_numTips], solutionDictionary) ; // We start obtaining keys at the root.
}

AugTree::AugTree(const umat & edgeMatrix, const vec & limProbs, const uint numTips, const uint rateCategIndex)
{
  _numTips = numTips ;
  _rateCateg = rateCategIndex ;
  _limProbs = limProbs ;
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
}

void AugTree::CopyAugTreeNonPointer(AugTree * sourceAugTree) 
{
  uint sourceVertexIndex = 0 ; // Defining an iterator on sourceAugTree->GetVertexVector() cannot be done without copying it, which I can't explain for now.
  for (auto & i : _vertexVector)
  {
    i->EnterCommonInfo(sourceAugTree->GetVertexVector().at(sourceVertexIndex)) ;
    i->EnterSolution(sourceAugTree->GetVertexVector().at(sourceVertexIndex)) ;
    sourceVertexIndex++ ;
  }
};

void AugTree::AssociateTransProbMatrices(const uvec & clusterMRCAs, const mat & withinTransProbMatrix, const mat & betweenTransProbMatrix) {
  
  // By default all nodes are considered between clusters.
  for (auto & i : _vertexVector)
  {
    i->SetTransProbMatrix(betweenTransProbMatrix, _rateCateg, false) ;
  }
  
  for (auto & i : clusterMRCAs)
  {
    if (i > _numTips) { // Again, clusterMRCAs is based on the R convention, hence >, and not >=.
      for (auto & j : _vertexVector[i-1]->GetChildren()) {
        BindMatrix(j, withinTransProbMatrix, true) ; 
      }  
    }
  }
}

void AugTree::BindMatrix(TreeNode * vertex, const mat & transProbMatrix, const bool withinCluster)
{
  vertex->SetTransProbMatrix(transProbMatrix, _rateCateg, withinCluster) ;
  if (!(vertex->GetChildren()[0] == NULL)) { // A null pointer indicates that we've reached an input node.
    for (auto & i : vertex->GetChildren())
    {
      BindMatrix(i, transProbMatrix, withinCluster) ;
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
  if (!vertex->CanFindKey()) // If a vertex can be solved, then its children have patterns assigned to them.
  {
    for (auto & i : vertex->GetChildren())
    {
      ComputeKeys(i, solutionDictionary) ;
    }
  }
  vertex->DeriveKey(solutionDictionary) ;
}

void AugTree::SolveRoot(solutionDictionaryType & solutionDictionary) {
  //PatternLookup(solutionDictionary, _vertexVector[_numTips]) ;
  TrySolve(_vertexVector[_numTips], solutionDictionary) ;
  _likelihood = dot(_vertexVector[_numTips]->GetSolution(), _limProbs) ;
}

void AugTree::InitializeVertices(const std::vector<uvec> & alignmentBinOneLocusOneRate)
{
  std::vector<TreeNode *>::iterator treeIter = _vertexVector.begin() ;
  for (auto & i : alignmentBinOneLocusOneRate)
  {
    (*treeIter)->SetInput(i) ; 
    (*treeIter)->ToggleSolved() ;
    treeIter++ ;
  }
}

void AugTree::PatternLookup(solutionDictionaryType & solutionDictionary, TreeNode * currentNode) {
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
      currentNode->ToggleSolved() ;
    }
  }
}

void AugTree::TrySolve(TreeNode * vertex, solutionDictionaryType & solutionDictionary)
{
  if (!(vertex->IsSolved())) // Could be solved because of the pattern lookup.
  {
    if (!vertex->CanSolve())
    {
      for (auto & i : vertex->GetChildren())
      {
        TrySolve(i, solutionDictionary) ;
      }
    }
    vertex->ComputeSolution(solutionDictionary) ;
  }
}

void AugTree::BindMatrixBetween(TreeNode * vertex, const mat & transProbMatrix)
{
  vertex->SetTransProbMatrix(transProbMatrix, _rateCateg, false) ;
  vertex->ToggleSolved() ;
  if (!(vertex->GetChildren().at(0) == NULL)) 
  {
    if (!(vertex->GetChildren().at(0)->GetWithinParentBranch())) 
    { // We're changing between-cluster transition probabilities only.
      for (auto & i : vertex->GetChildren()) 
      {
        BindMatrixBetween(i, transProbMatrix) ;
      }
    }
    else
    {
      vertex->ToggleSolved() ; // Not very clear... What happens here is that _solved is being set back to true for nodes supporting clusters. This is logical, since their solution is not assumed to have changed due to the change in the between-cluster transition probabilities.
    }
  }
}

void AugTree::InvalidateAll() // Assumes that the tree starts fully solved.
{
  for (auto & i : _vertexVector)
  {
    i->ToggleSolved() ;
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

Forest::Forest(const IntegerMatrix & edgeMatrix, const NumericVector & clusterMRCAs, const List & alignmentBin, const List & withinTransProbMatList, const List & betweenTransProbMatList, const NumericVector & limProbs, const uint numTips, const uint numLoci, solutionDictionaryType & solutionDictionary)
{
  _numLoci = numLoci ;
  _numRateCats = withinTransProbMatList.size() ;
  _solutionDictionary = solutionDictionary ;
  _forest.reserve(alignmentBin.size()*withinTransProbMatList.size()) ;
  std::vector<mat> withinTransProbsMats = as<std::vector<mat>>(withinTransProbMatList) ;
  std::vector<mat> betweenTransProbsMats = as<std::vector<mat>>(betweenTransProbMatList) ;
  umat edgeMatrixRecast = as<umat>(edgeMatrix) ;
  uvec clusterMRCAsRecast = as<uvec>(clusterMRCAs) ;
  Col<double> limProbsRecast = as<Col<double>>(limProbs) ;
  
  _randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ; // This is the random number generator. It's initialized when the Forest is built, and the seed is 0 by default.

  for (auto & i : alignmentBin) // Iterating on loci...
  {
    uint rateCategIndex = 0 ;
    uint locusRateIndex = 0 ;
    
    for (auto j : zip(withinTransProbsMats, betweenTransProbsMats))
    {
      AugTree * LocusRateAugTree = new AugTree(edgeMatrixRecast, clusterMRCAsRecast, j.get<0>(), j.get<1>(), as<std::vector<uvec>>(i), limProbsRecast, numTips, rateCategIndex, _solutionDictionary) ;
      _forest.push_back(LocusRateAugTree) ;
      rateCategIndex++ ;
    }
  }
}

Forest::Forest(const IntegerMatrix & edgeMatrix, const vec & limProbs, uint numRateCats, uint numLoci, uint numTips, gsl_rng * ranNumGenerator, solutionDictionaryType solutionDictionary)
{
  _numRateCats = numRateCats ;
  _numLoci = numLoci ;
  _randomNumGenerator = ranNumGenerator ;
  _forest.reserve(numLoci*numRateCats) ;
  umat edgeMatrixRecast = as<umat>(edgeMatrix) ;
  _solutionDictionary = solutionDictionary ;
  
  for (uint i = 0; i < numLoci; i++) // Iterating on loci...
  {
    uint rateCategIndex = 0 ;
    
    for (uint j = 0 ; j < numRateCats ; j++)
    {
      AugTree * LocusRateAugTree = new AugTree(edgeMatrixRecast, limProbs, numTips, j) ;
      _forest.push_back(LocusRateAugTree) ; 
      rateCategIndex++ ;
    }
  }
}


void Forest::ComputeLoglik()
{
  //#pragma omp parallel for 
  for (std::vector<AugTree *>::iterator forestIter = _forest.begin(); forestIter < _forest.end(); forestIter++) // This syntax is compatible with openMP, unlike the more conventional 'for (auto & i : myVec')
  {
    (*forestIter)->SolveRoot(_solutionDictionary) ;
  }
  
  // Now, we must average likelihoods across rate categories for each locus, log the output, and sum the resulting logs.
  Col<double> rateAveragedLogLiks(_numLoci) ;
  Col<double> likAcrossRatesLoci(_forest.size()) ;
  std::transform(_forest.begin(), _forest.end(), likAcrossRatesLoci.begin(), [] (AugTree * myTree) {return myTree->GetLikelihood() ;}) ;
  
  for (uint i = 0; i < rateAveragedLogLiks.size(); i++)
  {
    rateAveragedLogLiks[i] = log(mean(likAcrossRatesLoci.rows(_numRateCats*i, _numRateCats*(i+1) - 1))) ;
  }
  _loglik = sum(rateAveragedLogLiks) ;
}

void Forest::AmendBetweenTransProbs(std::vector<mat> & newBetweenTransProbs)
{
  uint rateCategIndex = 0 ;
  
  for (auto & tree : _forest)
  {
    tree->BindMatrixBetween(tree->GetVertexVector().at(tree->GetNumTips()), newBetweenTransProbs.at(rateCategIndex)) ;
    rateCategIndex = littleCycle(rateCategIndex+1, newBetweenTransProbs.size()) ;
  }
}

void Forest::AmendWithinTransProbs(std::vector<mat> & withinTransProbs, uvec & clusterMRCAs) 
{
  uint rateCategIndex = 0 ;
  for (auto &augtree : _forest)
  {
    for (auto & mrca : clusterMRCAs)
    {
      if (mrca <= augtree->GetNumTips())
      {
        augtree->BindMatrix(augtree->GetVertexVector().at(mrca - 1), withinTransProbs.at(rateCategIndex), true) ;
      }
    } 
    augtree->InvalidateAll() ;
    rateCategIndex = littleCycle(rateCategIndex, withinTransProbs.size()) ;
  }
}

void Forest::HandleSplit(uint clusMRCAtoSplit, std::vector<mat> & betweenTransProbsMats)
{
  uint rateCateg = 0 ;
  
  for (auto & augtree : _forest)
  {
    augtree->GetVertexVector().at(clusMRCAtoSplit - 1)->InvalidateSolution() ;
    for (auto & childNode : augtree->GetVertexVector().at(clusMRCAtoSplit - 1)->GetChildren())
    {
      childNode->SetTransProbMatrix(betweenTransProbsMats.at(rateCateg), rateCateg, false) ;
    }
    rateCateg = littleCycle(rateCateg + 1, betweenTransProbsMats.size()) ;
  }
}

void Forest::HandleMerge(uvec & clusMRCAstoMerge, std::vector<mat> & withinTransProbsMats)
{
  uint rateCateg = 0 ;
  for (auto & augtree : _forest)
  {
    augtree->GetVertexVector().at(clusMRCAstoMerge.at(0) - 1)->GetParent()->InvalidateSolution() ; // Elements of clusMRCAsToMerge should all have the same parent to allow a merge to occur.
    for (auto & oldClusterMRCA : clusMRCAstoMerge)
    {
      augtree->GetVertexVector().at(oldClusterMRCA - 1)->SetTransProbMatrix(withinTransProbsMats.at(rateCateg), rateCateg, true) ;
    }
    rateCateg = littleCycle(rateCateg + 1, withinTransProbsMats.size()) ;
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
  uint originAugTreeIndex = 0 ;
  for (auto & i : _forest)  {
    
    i->CopyAugTreeNonPointer(originForest->GetForest().at(originAugTreeIndex)) ;
    originAugTreeIndex++ ;
  }
}