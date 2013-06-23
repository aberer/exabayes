#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "LnlRestorer.hpp"

class LikelihoodEvaluator
{
public: 
  LikelihoodEvaluator(LnlRestorerPtr restorer); 
  
  double evaluatePartitions( TreeAln &traln, const vector<nat>& partitions)  ; 
  void evalSubtree( TreeAln &traln, const Branch &evalBranch)    ; 
  double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 
  void findVirtualRoot(const TreeAln &traln, Branch &result) const ; 

  void imprint(const TreeAln &traln) { restorer->resetRestorer(traln); }
  void resetToImprinted(TreeAln &traln){restorer->restoreArrays(traln); }

  static bool disorientNode( nodeptr p); 
  static void disorientTree(TreeAln &traln, const Branch &root) ; 
  static void disorientSubtree(TreeAln &traln, const Branch &branch) ; 


private: 			// METHODS
  void coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked); 



private: 			// ATTRIBUTES
  LnlRestorerPtr restorer;    

}; 

typedef shared_ptr<LikelihoodEvaluator> LikelihoodEvaluatorPtr; 

#endif

