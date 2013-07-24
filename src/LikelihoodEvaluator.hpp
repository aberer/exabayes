#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include "axml.h"
#include "TreeAln.hpp"
#include "Branch.hpp"
#include "LnlRestorer.hpp"

class LikelihoodEvaluator
{
public: 
  LikelihoodEvaluator(){}
  virtual ~LikelihoodEvaluator() {} 

  /**
     @brief: evaluate a list of partitions. This is always a full traversal 
   */
  virtual double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions, bool fullTraversal)  = 0; 
  /** 
      @brief evaluate a subtree. This used to be the
      newview-command. The this-node is the important part, the
      that-node of the branch only specifies the orientation.
   */ 
  virtual void evalSubtree( TreeAln &traln, const Branch &evalBranch)  = 0   ; 
  /** 
      @brief evaluation at a given branch  
   */ 
  virtual double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal ) = 0 ; 
  /** 
      @brief find the root branch in the current tree    
   */
  Branch findVirtualRoot(const TreeAln &traln) const ;   
  /** 
      @brief make the current state in the tree resettable (only needed for chain.cpp)
   */ 
  virtual void imprint(const TreeAln &traln) = 0; 
  /**
     @brief use backup likelihood arrays to reset the state 
   */ 
  virtual void resetToImprinted(TreeAln &traln) = 0;  
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

  virtual std::unique_ptr<LikelihoodEvaluator> clone() const = 0; 

  void expensiveVerify(TreeAln &traln);
  void setDebugTraln(std::shared_ptr<TreeAln> _debugTraln); 

protected: 			// METHODS
  void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal); 
  void coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked); 

protected: 			// ATTRIBUTES
  std::shared_ptr<TreeAln> debugTraln;  
  bool verifyLnl; 
  double prevLnl; 
  std::vector<double> partitionLnls;  
}; 

#endif

