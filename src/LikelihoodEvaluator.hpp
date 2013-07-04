#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "LnlRestorer.hpp"

class LikelihoodEvaluator
{
public: 
  LikelihoodEvaluator(std::shared_ptr<LnlRestorer> restorer); 

  double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions)  ; 
  void evalSubtree( TreeAln &traln, const Branch &evalBranch)    ; 
  double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 
  void findVirtualRoot(const TreeAln &traln, Branch &result) const ; 

  void imprint(const TreeAln &traln) { restorer->resetRestorer(traln); }
  void resetToImprinted(TreeAln &traln){restorer->restoreArrays(traln); }

  static bool disorientNode( nodeptr p); 
  static void disorientTree(TreeAln &traln, const Branch &root) ; 
  static void disorientSubtree(TreeAln &traln, const Branch &branch) ; 

  void evaluateFullNoBackup(TreeAln& traln); 
  void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal); 

#ifdef DEBUG_LNL_VERIFY
  // BAD
  void expensiveVerify(TreeAln &traln);
  void setDebugTraln(std::shared_ptr<TreeAln> _debugTraln); 
#endif

private: 			// METHODS
  void coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked); 


private: 			// ATTRIBUTES
  std::shared_ptr<LnlRestorer> restorer;    
#ifdef DEBUG_LNL_VERIFY
  std::shared_ptr<TreeAln> debugTraln;  
  bool verifyLnl; 
#endif

}; 




#endif

