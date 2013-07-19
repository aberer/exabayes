#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "LnlRestorer.hpp"

class LikelihoodEvaluator
{
public: 
  /**
     @brief: evaluate a list of partitions. This is always a full traversal 
   */
  virtual double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions) = 0; 
  /** 
      @brief evaluate a subtree. This used to be the
      newview-command. The this-node is the important part, the
      that-node of the branch only specifies the orientation.
   */ 
  virtual void evalSubtree( TreeAln &traln, const Branch &evalBranch)  = 0   ; 
  /** 
      @brief evaluation at a given branch  
   */ 
  virtual double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal ) = 0; 
  /** 
      @brief find the root branch in the current tree    
   */
  Branch findVirtualRoot(const TreeAln &traln) const ;   
  /** 
      @brief make the current state in the tree resettable (only needed for chain.cpp)
   */ 
  void imprint(const TreeAln &traln); 
  /**
     @brief use backup likelihood arrays to reset the state 
   */ 
  void resetToImprinted(TreeAln &traln) ; 
  /** 
      @brief invalidate the orientation at a given node 
   */ 
  static bool disorientNode( nodeptr p); 
  /**
     @brief destroy the orientation of the entire tree (e.g., to enforce re-evaluation)
   */
  static void disorientTree(TreeAln &traln, const Branch &root) ; 
  /** 
      @brief destroy the orientation of a subtree 
   */ 
  static void disorientSubtree(TreeAln &traln, const Branch &branch) ; 


  /** 
      @brief conduct full traversal and evaluation on tree, do not use
      likelihood backup arrays (rarely useful)
   */ 
  void evaluateFullNoBackup(TreeAln& traln); 

#ifdef DEBUG_LNL_VERIFY
  // BAD
  void expensiveVerify(TreeAln &traln);
  void setDebugTraln(std::shared_ptr<TreeAln> _debugTraln); 
#endif

protected: 			// METHODS
  void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal); 
  void coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked); 

protected: 			// ATTRIBUTES
  std::shared_ptr<LnlRestorer> restorer;    
#ifdef DEBUG_LNL_VERIFY
  std::shared_ptr<TreeAln> debugTraln;  
  bool verifyLnl; 
#endif

  double prevLnl; 
}; 

#endif

