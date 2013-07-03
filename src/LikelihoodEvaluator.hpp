#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "LnlRestorer.hpp"

class LikelihoodEvaluator
{
public: 
  LikelihoodEvaluator(shared_ptr<LnlRestorer> restorer); 

  double evaluatePartitions( TreeAln &traln, const vector<nat>& partitions)  ; 
  void evalSubtree( TreeAln &traln, const Branch &evalBranch)    ; 
  double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 
  void findVirtualRoot(const TreeAln &traln, Branch &result) const ; 

  void imprint(const TreeAln &traln) { restorer->resetRestorer(traln); }
  void resetToImprinted(TreeAln &traln){restorer->restoreArrays(traln); }

  static bool disorientNode( nodeptr p); 
  static void disorientTree(TreeAln &traln, const Branch &root) ; 
  static void disorientSubtree(TreeAln &traln, const Branch &branch) ; 

  // BAD
  void evaluateFullNoBackup(TreeAln& traln)
  {
#ifdef DEBUG_EVAL  
    cout << "conducting full evaluation, no backup created" << endl; 
#endif
  
    exa_evaluateGeneric(traln,traln.getTr()->start,TRUE );   
    expensiveVerify(traln);
  }

void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL != 0
  evaluateGeneric(traln.getTr(), traln.getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(traln.getTr(), start, fullTraversal); 
#endif  
}
  // BAD
  void expensiveVerify(TreeAln &traln);

private: 			// METHODS
  void coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked); 


private: 			// ATTRIBUTES
  shared_ptr<LnlRestorer> restorer;    

}; 


#endif

