#include "AugTree.h"
#include "InputNode.h"
#include "IntermediateNode.h"
#include "TreeNode.h"
#include <unordered_map>
#include "helper.h"

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(const umat & edgeMatrix, const uvec & clusterMRCAs, std::vector<uvec> * alignmentBinOneLocusOneRate, const uint rateCategIndex, solutionDictionaryType & solutionDictionary, const uint & numTips)
{ 
  _exponentContainer = 0 ;
  _rateCateg = rateCategIndex ;
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy, numTips) ;
  InitializeVertices(alignmentBinOneLocusOneRate, solutionDictionary) ;
  AssociateTransProbMatrices(clusterMRCAs, numTips) ;
}

AugTree::AugTree(const umat & edgeMatrix, const uint & rateCategIndex, const uint & numTips)
{
  _exponentContainer = 0 ;
  _rateCateg = rateCategIndex ;
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy, numTips) ;
}

void AugTree::CopyAugTreeNonPointer(AugTree * sourceAugTree) 
{
  _exponentContainer = sourceAugTree->GetExponentContainer() ;
  uint sourceVertexIndex = 0 ; // Defining an iterator on sourceAugTree->GetVertexVector() cannot be done without copying it, which I can't explain for now.
  for (auto & i : _vertexVector)
  {
    i->EnterCommonInfo(sourceAugTree->GetVertexVector().at(sourceVertexIndex)) ;
    i->EnterInput(sourceAugTree->GetVertexVector().at(sourceVertexIndex)) ;
    sourceVertexIndex++ ;
  }
}

void AugTree::AssociateTransProbMatrices(const uvec & clusterMRCAs, const uint & numTips) 
{
  // By default all nodes are considered between clusters.
  for (auto & i : _vertexVector)
  {
    i->SetWithinParentBranch(false) ;
  }
  
  for (auto & i : clusterMRCAs)
  {
    if (i > numTips) { // Again, clusterMRCAs is based on the R convention, hence >, and not >=.
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

void AugTree::BuildTree(const umat & edgeMatrix, const uint & numTips)
{
  _vertexVector.reserve(edgeMatrix.n_rows + 1) ;

  // We create the tips. Note that tip 1 should correspond to vertex 1 in the original (the one in the phylo object) edgeMatrix

  for (uint i = 0; i < numTips; i++) {
    InputNode * newNode = new InputNode{} ;
    _vertexVector.push_back(newNode) ;
  } ;

  // We add the internal nodes
  for (uint i = 0 ; i < edgeMatrix.n_rows - numTips + 1; i++) {
    IntermediateNode * newNode = new IntermediateNode{};
    _vertexVector.push_back(newNode) ; 
  } ;
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _vertexVector.size(); i++) {
    _vertexVector[i]->SetId(i) ; 
  } ;
  // The vertices are all disjoint, the next loop defines their relationships
  // The iterator follows columns.
  for (umat::const_iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end()-edgeMatrix.n_rows; iter++)
  {
    _vertexVector[*iter]->AddChild(_vertexVector[*(iter+edgeMatrix.n_rows)]) ;
    _vertexVector[*(iter+edgeMatrix.n_rows)]->SetParent(_vertexVector[*iter]) ;
  }
}

void AugTree::BuildTreeNoAssign(const umat & edgeMatrix)
{
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _vertexVector.size(); i++) {
    _vertexVector[i]->RemoveChildren() ; // Might be more memory-efficient to overwrite the children... 
  } ;
  // The vertices are all disjoint, the next loop defines their relationships
  // The iterator follows columns.
  for (umat::const_iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end()-edgeMatrix.n_rows; iter++)
  {
    _vertexVector[*iter]->AddChild(_vertexVector[*(iter+edgeMatrix.n_rows)]) ;
    _vertexVector[*(iter+edgeMatrix.n_rows)]->SetParent(_vertexVector[*iter]) ;
  }
}



void AugTree::ComputeKeys(TreeNode * vertex, solutionDictionaryType & solutionDictionary, const uint withinMatListIndex, const uint betweenMatListIndex)
{
  if (!vertex->IsKeyDefined())
  {
    if (!vertex->CanFindKey()) // If a vertex can be solved, then its children have patterns assigned to them.
    {
      for (auto & i : vertex->GetChildren())
      {
        ComputeKeys(i, solutionDictionary, withinMatListIndex, betweenMatListIndex) ;
      }
    }
    uint matListIndex = betweenMatListIndex ; 
    if (vertex->GetChildren().at(0)->GetWithinParentBranch())
    {
      matListIndex = withinMatListIndex ;
    }
    vertex->DeriveKey(solutionDictionary, _rateCateg, matListIndex) ;
  }
}

void AugTree::SolveRoot(solutionDictionaryType & solutionDictionary, const mat & withinTransProbMat, const mat & betweenTransProbMat, const vec & limProbs, const uint & numTips)
{
  TrySolve(_vertexVector[numTips], solutionDictionary, withinTransProbMat, betweenTransProbMat) ;
  _likelihoodProp = dot(_vertexVector[numTips]->GetSolution(solutionDictionary, _rateCateg), limProbs) ;
}

void AugTree::InitializeVertices(std::vector<uvec> * alignmentBinOneLocusOneRate, solutionDictionaryType & solutionDictionary)
{ 
  std::vector<TreeNode *>::iterator treeIter = _vertexVector.begin() ;
  for (auto & i : *alignmentBinOneLocusOneRate)
  {
    (*treeIter)->SetInput(&i) ;
    (*treeIter)->DeriveKey(solutionDictionary, _rateCateg, 0) ; // The 0 is a placeholder. It is not used in deriving the key for the tip configurations. 
    treeIter++ ;
  }
}

void AugTree::PatternLookup(solutionDictionaryType & solutionDictionary, TreeNode * currentNode) 
{
  if (!currentNode->IsSolved() && !solutionDictionary->empty())
  { // If the node is already solved, no need to update it with a stored pattern.
    if (solutionDictionary->at(_rateCateg).count(currentNode->GetDictionaryKey()) == 0)
    {
      for (auto & i : currentNode->GetChildren())
      {
        PatternLookup(solutionDictionary, i) ;
      }
    }
    else
    {
      //currentNode->SetSolution((*solutionDictionary)[currentNode->GetDictionaryKey()]) ;
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
      vertex->ComputeSolution(solutionDictionary, withinTransProbMat, &_exponentContainer, _rateCateg) ;
    }
    else
    {
      vertex->ComputeSolution(solutionDictionary, betweenTransProbMat, &_exponentContainer, _rateCateg) ;
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

void AugTree::RearrangeTreeNNI(uint vertexId1, uint vertexId2, solutionDictionaryType solutionDictionary) 
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

umat AugTree::BuildEdgeMatrix(const uint & numTips)
{
  umat edgeMatrix(_vertexVector.size()-1, 2, fill::zeros) ;
  uint lineNum = 0 ;
  AddEdgeRecursion(edgeMatrix, lineNum, _vertexVector.at(numTips)) ;
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
