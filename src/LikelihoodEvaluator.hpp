#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"

class LikelihoodEvaluator
{
public: 
  double evaluatePartitions( TreeAln &traln, const vector<nat>& partitions)  ; 
  void evalSubtree( TreeAln &traln, const Branch &evalBranch)    ; 
  double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 

  void findVirtualRoot(const TreeAln &traln, Branch &result) const ; 

  static bool disorientNode( nodeptr p); 
  static void disorientTree(TreeAln &traln, const Branch &root) ; 
  static void disorientSubtree(TreeAln &traln, const Branch &branch) ; 

private: 
  // TODO lnl restorer 

}; 



#endif

