#include "AugTree.h"
#include "IntermediateNode.h"
#include "InputNode.h"
#include<vector.h>
#include<memory.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

using namespace Rcpp ;
using namespace arma ;

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}


AugTree::AugTree(IntegerMatrix & edgeMatrix, List & alignmentBinOneLocus, List & withinTransProbMatList, List & betweenTransProbMatList, NumericVector & limProbs) 
{
  _limProbs = as<vec>(limProbs) ;
  _withinTransProbMatrixVec = as<std::vector<mat>>(withinTransProbMatList) ;
  _betweenTransProbMatrixVec = as<std::vector<mat>>(betweenTransProbMatList) ;
  _alignmentBin = as<std::vector<longVec>>(alignmentBin) ;
  umat edgeMatrixCopy(as<umat>(edgeMatrix)) ;
  edgeMatrixCopy = edgeMatrixCopy - 1 ;
  _numRateCats = _withinTransProbMatrixVec.size() ;
  BuildTree(edgeMatrixCopy) ;
  SetPatterns() ;
}

void AugTree::BuildTree(umat & edgeMatrix) 
{
  
  _tree.reserve(edgeMatrix.n_rows + 1) ;
  
  // We create the tips. Note that tip 1 should correspond to vertex 1 in the original (the one in the phylo object) edgeMatrix
  for (std::vector<longVec>::iterator iter = _alignmentBin.begin(); iter < _alignmentBin.end() ; iter++) {
    InputNode newNode(*iter, _numRateCats) ;  
    _tree.push_back(&newNode) ;
  } ;
  
  // We add the internal nodes
  for (uint i = 0 ; i < edgeMatrix.n_rows - _alignmentBin.size() + 1; i++) {
    IntermediateNode myNode{} ;
    _tree.at(i) = &myNode ;
  } ;
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _tree.size(); i++) {
    _tree.at(i)->SetId(i+1) ; // R want edge labels that start at 1.
  } ;
  
  // The vertices are all disjoint, the next loop defines their relationships
  for (umat::iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end(); iter = iter + 2) 
  { 
    _tree.at(*iter)->AddChild(_tree.at(*(iter+1))) ;
    _tree.at(*(iter+1))->SetParent(_tree.at(*iter)) ;
  }
}

SEXP AugTree::BuildEdgeMatrix() 
{
  //TO_DO
}

void AugTree::InitializePatterns(TreeNode * startNode) 
{
  for (auto & i : startNode->_children)
  {
    if (i->CanSolve())
    {
      i->SetPattern() ;
      i->ToggleSolved() ;
    }
    else
    {
     InitializePatterns(i) ;
    }
  }
  // The nodes had been marked as solved only so that we could get a pattern indicator for each. They are not really solved and so, must be marked as such.
  // Note that, as it should, ToggleSolved() does nothing for InputNodes.
  for (auto & i : _tree)
  {
    i->ToggleSolved() ;
  }
}

bool AugTree::SolveRoot() {
  //TreeNode * currentPoint = _tree.at(_alignmentBin.size()) ; // _alignmentBin has 2 levels: the outer level is node, the inner level is locus.
  for (uint locusNum = 0; locusNum < _alignmentBin[0].size(); locusNum++) 
  {
    for (uint rateNum = 0; rateNum < _numRateCats; rateNum++) 
    {
      InitializeFromDictionary(rateNum, locusNum) ;
      SolveOneRateOneLocus(rateNum, locusNum) ; 
    }
    
  }
}

void setPatterns {
  //TO_DO
}

void AugTree::InitializeFromDictionary(uint rateNum, uint locusNum) 
{
  TreeNode * currentPoint = _tree.at(_alignmentBin.size()) ;
  if (_dictionary.find(currentPoint->GetPattern().at(rateNum).at(locusNum)) != _dictionary.end()) 
  {
    
  }
  
        
  
}