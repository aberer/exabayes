#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include <vector>

#include "axml.h"
#include "TreeAln.hpp"
#include "ArrayRestorer.hpp"

class LikelihoodEvaluator
{
public: 
  LikelihoodEvaluator(){}
  virtual ~LikelihoodEvaluator() {} 

  /**
     @brief: evaluate a list of partitions. This is always a full traversal 
   */
  double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions, bool fullTraversal)  ;
  virtual double evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain& root,  const std::vector<nat>& partitions, bool fullTraversal)  = 0; 
  /** 
      @brief evaluate a subtree. This used to be the
      newview-command. The this-node is the important part, the
      that-node of the branch only specifies the orientation.
   */ 
  virtual void evalSubtree( TreeAln &traln, const BranchPlain &evalBranch)  = 0   ; 
  /** 
      @brief evaluation at a given branch  
   */ 
  virtual double evaluate(TreeAln &traln, const BranchPlain &evalBranch,  bool fullTraversal ) = 0 ; 
  /** 
      @brief find the root branch in the current tree    
   */
  BranchPlain findVirtualRoot(const TreeAln &traln) const ;   
  /** 
      @brief mark a node as dirty  
   */ 
  void markDirty(nat partition, nat nodeId); 
  /** 
      @brief marks a node as clean
   */ 
  void markClean(nat partition, nat nodeId); 

  /** 
      @brief make the current state in the tree resettable (only needed for chain.cpp)
   */ 
  virtual void imprint(const TreeAln &traln) = 0; 
  /**
     @brief use backup likelihood arrays to reset the state 
   */ 
  virtual void resetToImprinted(TreeAln &traln) = 0;    
  /** 
      @brief use backup likelihood
   */ 
  virtual void resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions) = 0; 
  /** 
      @brief invalidate the orientation at a given node 
   */ 
  static bool disorientNode( nodeptr p); 
  /**
     @brief destroy the orientation of the entire tree (e.g., to enforce re-evaluation)
   */
  static void disorientTree(TreeAln &traln, const BranchPlain &root) ; 
  /** 
      @brief destroy the orientation of a subtree 
   */ 
  static void disorientSubtree(TreeAln &traln, const BranchPlain &branch) ; 

  virtual std::unique_ptr<LikelihoodEvaluator> clone() const = 0; 

  void expensiveVerify(TreeAln &traln);
  void setDebugTraln(std::shared_ptr<TreeAln> _debugTraln); 
  /** 
      @brief sets a flag that a likelihood array for a node and a
      partition is invalid
  */ 
  void setDirty(nat partition, nat nodeId, bool isDirty); 
  /** 
      @brief indicates whether a likelihood array for a given node and
      partition is invalid
  */ 
  bool isDirty(nat partition,nat nodeId)  const ; 

protected: 			// METHODS
  void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal); 
  void coreEvalSubTree(TreeAln& traln, nodeptr p, boolean masked); 

protected: 			// ATTRIBUTES
  std::shared_ptr<TreeAln> debugTraln;  
  bool verifyLnl; 
  double prevLnl; 
  std::vector<double> partitionLnls;  
  std::vector<std::vector<bool> > dirty; 
}; 

#endif

