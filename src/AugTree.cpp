#include "AugTree.h"
#include "InputNode.h"
#include "IntermediateNode.h"
#include "TreeNode.h"
#include <unordered_map>
#include "helper.h"

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(const umat & edgeMatrix, const uvec & clusterMRCAs, std::vector<std::vector<uvec>> * alignmentBin, solutionDictionaryType & solutionDictionary, const uint & withinMatListIndex, const uint & betweenMatListIndex, const uint & numRateCats, gsl_rng * RNGpoint, unsigned int & numThreads, boost::asio::io_service * ioService, boost::asio::io_service::work * myWork) : _workObject(myWork), _ioService(ioService)
{
  _solutionDictionary = solutionDictionary ;
  _alignmentBinReference = alignmentBin ;
  _numTips = alignmentBin->size() ;
  _numLoci = alignmentBin->at(0).size() ;
  _numRateCats = numRateCats ;
  _numThreads = numThreads ;
  
  _withinMatListIndex = withinMatListIndex ;
  _betweenMatListIndex = betweenMatListIndex ;
  _randomNumGenerator = RNGpoint ;
  
  umat edgeMatrixCopy(edgeMatrix) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  BuildTree(edgeMatrixCopy) ;
  InitializeVertices() ;
  AssociateTransProbMatrices(clusterMRCAs) ;
  
  /*
   * This will add numThreads threads to the thread pool. I wonder if the master thread should also perform some work...
   */
  for (uint i = 0; i < numThreads; i++)
  {
    _threadpool.create_thread(
    boost::bind(&boost::asio::io_service::run, _ioService)
  );
  }
}

void AugTree::AssociateTransProbMatrices(const uvec & clusterMRCAs) 
{
  // By default all nodes are considered between clusters.
  for (auto & i : _vertexVector)
  {
    i->SetWithinParentBranch(false) ;
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

void AugTree::BuildTree(const umat & edgeMatrix)
{
  _vertexVector.reserve(edgeMatrix.n_rows + 1) ;

  // We create the tips. Note that tip 1 should correspond to vertex 1 in the original (the one in the phylo object) edgeMatrix

  for (uint i = 0; i < _numTips; i++) {
    InputNode * newNode = new InputNode(_numLoci, _numRateCats, _solutionDictionary) ;
    _vertexVector.push_back(newNode) ;
  }

  // We add the internal nodes
  for (uint i = 0 ; i < edgeMatrix.n_rows - _numTips + 1; i++) {
    IntermediateNode * newNode = new IntermediateNode(_numLoci, _numRateCats, _solutionDictionary) ;
    _vertexVector.push_back(newNode) ; 
  }
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _vertexVector.size(); i++) {
    _vertexVector[i]->SetId(i) ; 
  }
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
    _vertexVector.at(i)->RemoveChildren() ; // Might be more memory-efficient to overwrite the children... 
  } ;
  // The vertices are all disjoint, the next loop defines their relationships
  // The iterator follows columns.
  for (umat::const_iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end()-edgeMatrix.n_rows; iter++)
  {
    _vertexVector.at(*iter)->AddChild(_vertexVector.at(*(iter+edgeMatrix.n_rows))) ;
    _vertexVector.at(*(iter+edgeMatrix.n_rows))->SetParent(_vertexVector.at(*iter)) ;
  }
}

void AugTree::InitializeVertices()
{ 
  std::vector<TreeNode *>::iterator treeIter = _vertexVector.begin() ;
  for (auto & i : *_alignmentBinReference)
  {
    (*treeIter)->SetInput(&i) ;
    (*treeIter)->InitMapAndIterVec(_solutionDictionary, _numRateCats) ; // The map will have containers for the input nodes, but since their solutions are known from the start, we just need to put them in the map verbatim.
    treeIter++ ;
  }
}

void AugTree::TrySolve(TreeNode * vertex, const std::vector<mat> & withinTransProbMats, const std::vector<mat> & betweenTransProbMats, const uint & locusNum)
{
  if (!(vertex->IsSolved()))
  {
    if (!vertex->CanSolve())
    {
      for (auto & i : vertex->GetChildren())
      {
        TrySolve(i, withinTransProbMats, betweenTransProbMats, locusNum) ;
      }
    }
    if (vertex->GetChildren().at(0)->GetWithinParentBranch()) 
    {  // This junction is within a cluster. 
      vertex->ComputeSolution(_solutionDictionary, withinTransProbMats, locusNum, _withinMatListIndex, _mutex) ;
    }
    else
    {
      vertex->ComputeSolution(_solutionDictionary, betweenTransProbMats, locusNum, _betweenMatListIndex, _mutex) ;
    }
    if (locusNum == _numLoci-1)
    {
      vertex->SetSolved(true) ;
      vertex->SetUpdate(true) ;
    }
  }
}

void AugTree::InvalidateAll() // Assumes that the tree starts fully solved.
{
  for (auto & i : _vertexVector)
  {
    i->SetSolved(false) ;
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

std::vector<uint> AugTree::GetTwoVerticesForNNI(TreeNode * subtreeRoot, uvec & clusterMRCAs)
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
      selectedChildIndex = gsl_rng_uniform_int(_randomNumGenerator, i->GetChildren().size()) ;
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
      
      for (auto & child : currentVertex->GetChildren())
      {
        CheckAndInvalidateBetweenRecursive(child) ;
      }
    }
  }
}

void AugTree::ComputeLoglik(const std::vector<mat> & withinClusTransProbs, const std::vector<mat> & betweenClusTransProbs, const vec & limProbs)
{
  uint numElements = _numLoci*_numRateCats ;
  
  for (uint locusNum = 0 ; locusNum < _numLoci ; locusNum++) {
    TrySolve(_vertexVector[_numTips], withinClusTransProbs, betweenClusTransProbs, locusNum) ;
  }
  vec likPropVec(numElements, fill::zeros) ;
  
  for (uint locusIndex = 0 ; locusIndex < _numLoci ; locusIndex++) 
  {
    for (uint rateIndex = 0 ; rateIndex < _numRateCats ; rateIndex++)
    {
      likPropVec.at(locusIndex*_numRateCats + rateIndex) = dot(_vertexVector.at(_numTips)->GetSolution(locusIndex, rateIndex, _mutex), limProbs) ; 
    }
  }
  // Now, we must average likelihoods across rate categories for each locus, log the output, and sum the resulting logs.
  vec rateAveragedLogLiks(_numLoci, fill::zeros) ;
  fvec exponentVec(_numRateCats, fill::zeros) ;
  for (uint locusIndex = 0; locusIndex < _numLoci; locusIndex++)
  {
    for (uint rateIndex = 0; rateIndex < _numRateCats; rateIndex++)
    {
      for (uint i = _numTips ; i < _vertexVector.size() ; i++) // The first _numTips vertices are input nodes and therefore, their exponent factorization is trivially 0.
      {
        exponentVec.at(rateIndex)+=_vertexVector.at(i)->GetExponent(locusIndex, rateIndex) ;
      }
    }
    
    double maxExponent = max(exponentVec) ;
    exponentVec-=maxExponent ;
    
    double vectorIndex = _numRateCats*locusIndex ;
    double endVectorIndex = _numRateCats*(locusIndex+1) - 1 ;
    rateAveragedLogLiks.at(locusIndex) = log(mean(likPropVec.rows(vectorIndex, endVectorIndex)%exp(exponentVec))) + maxExponent;
  }
  
  _logLik = sum(rateAveragedLogLiks) ;
}

void AugTree::HandleSplit(uint clusMRCAtoSplit)
{
  _vertexVector.at(clusMRCAtoSplit - 1)->InvalidateSolution() ;
  
  for (auto & childNode : _vertexVector.at(clusMRCAtoSplit - 1)->GetChildren())
  {
    childNode->SetWithinParentBranch(false) ;
  }
}

void AugTree::HandleMerge(uvec & clusMRCAstoMerge)
{
  _vertexVector.at(clusMRCAstoMerge.at(0) - 1)->GetParent()->InvalidateSolution() ; // Elements of clusMRCAsToMerge should all have the same parent to allow a merge to occur.
  
  for (auto & oldClusterMRCA : clusMRCAstoMerge)
  {
    _vertexVector.at(oldClusterMRCA - 1)->SetWithinParentBranch(true) ;
  }
}

void AugTree::RestorePreviousConfig(const IntegerMatrix & edgeMat, const bool NNImoveFlag, const int & withinMatListIndex, const int & betweenMatListIndex, const IntegerVector & clusterMRCAs, bool splitMergeMove) 
{
  _withinMatListIndex = withinMatListIndex ;
  _betweenMatListIndex = betweenMatListIndex ;
 
  if (NNImoveFlag)
  {
    BuildTreeNoAssign(as<umat>(edgeMat)-1) ;
  }
  
  for (auto & vertex : _vertexVector)
  {
    vertex->RestoreIterVecAndExp() ;
  }
  if (splitMergeMove)
  {
    AssociateTransProbMatrices(as<uvec>(clusterMRCAs)) ; 
  }
}

void AugTree::NegateAllUpdateFlags()
{
  for (auto & i : _vertexVector)
  {
    i->NegateFlag() ;
  }
}

void AugTree::PrintSolutions(const uint & elementNum)
{
  // for (auto & i : _vertexVector)
  // {
  //   cout << "This is node " << i->GetId() << ".";
  //   cout << "Is my supporting branch within-cluster? " << i->GetWithinParentBranch() << ".\n";
  //   i->GetDictionaryIterator(elementNum, _numRateCats)->second.first.print("Solution from dictionary (could be rescaled):") ;
  // }
}