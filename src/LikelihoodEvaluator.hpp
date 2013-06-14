#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"

class LikelihoodEvaluator
{
public: 
  double evaluatePartitions( TreeAln &traln, const vector<nat>& partitions, bool fullTraversal)  ; 
  void evalSubtree( TreeAln &traln, const Branch &evalBranch)    ; 
  double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal )  ; 

  void findVirtualRoot(const TreeAln &traln, Branch &result) const ; 

  void disorientTree(TreeAln &traln, const Branch &root) const; 
  void disorientSubtree(TreeAln &traln, const Branch &branch) const; 

private: 
  // TODO lnl restorer 

}; 



#endif

