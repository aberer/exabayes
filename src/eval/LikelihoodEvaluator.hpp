#ifndef _TREE_EVALUATOR_HPP
#define _TREE_EVALUATOR_HPP

#include <vector>

#include "ArrayPolicy.hpp"
#include "axml.h"
#include "TreeAln.hpp"
#include "ArrayOrientation.hpp"


class LikelihoodEvaluator
{
public: 
  LikelihoodEvaluator(const TreeAln &traln, ArrayPolicy *plcy, std::shared_ptr<ArrayReservoir> arrayReservoir); 
  LikelihoodEvaluator( LikelihoodEvaluator &&rhs); 
  LikelihoodEvaluator( const LikelihoodEvaluator &rhs ); 
  LikelihoodEvaluator& operator=(LikelihoodEvaluator rhs) ; 
  friend void swap(LikelihoodEvaluator &lhs, LikelihoodEvaluator  &rhs); 

  /**
     @brief: evaluate a list of partitions. This is always a full traversal 
   */
  void evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain& root,  const std::vector<nat>& partitions, bool fullTraversal, bool debug)  ; 
  /** 
      @brief evaluate a subtree. This used to be the
      newview-command. The this-node is the important part, the
      that-node of the branch only specifies the orientation.
   */ 
  void evalSubtree(TreeAln  &traln, nat partition, const BranchPlain &evalBranch) ; 
  /** 
      @brief evaluation at a given branch  
   */ 
  void evaluate(TreeAln &traln, const BranchPlain &evalBranch,  bool fullTraversal, bool debug )  ; 
  /** 
      @brief mark a node as dirty  
   */ 
  void markDirty(const TreeAln &traln, nat nodeId); 
  void markDirty(const TreeAln &traln, nat partitionId, nat nodeId) ; 
  
  void debugPrintToCompute(const TreeAln &traln, const BranchPlain &root); 
  void debugPrintToComputeHelper(const TreeAln &traln, const BranchPlain &root ); 

  /** 
      @brief marks the entire tree dirty 
   */ 
  void markPartitionDirty(const TreeAln &traln, nat partition ); 
  /** 
      @brief indicates whether a likelihood array for a given node and
      partition is invalid
  */ 
  bool isDirty(nat partition,nat nodeId)  const ; 
  /** 
      @brief make the current state in the tree resettable (only needed for chain.cpp)
   */ 
  void imprint(const TreeAln &traln); 
  /** 
      @brief the main point of this is to free remaining memory  
   */ 
  void freeMemory() ; 
  /** 
      @brief use backup likelihood
   */ 
  void accountForRejection(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat> &invalidNodes); 

  ArrayOrientation getOrientation() const {return arrayOrientation;   }

  void expensiveVerify(TreeAln &traln,  BranchPlain root, double toVerify);
  void setDebugTraln(std::shared_ptr<TreeAln> _debugTraln); 

  ArrayReservoir& getArrayReservoir() {return *(_arrayReservoir) ; }

private: 			// METHODS 
  /** 
      @brief traverses the tree and reorients nodes to enforce
      re-computation
      
      @param the new virtual root
      
      @notice an important assumption is that the new virtual root is
      adjacent to invalidated likelihood arrays. Everything else would
      be extremely inefficient.
   */ 
  bool applyDirtynessToSubtree(TreeAln &traln, nat partition, const BranchPlain &branch); 
  void disorientDebug(TreeAln &traln, const BranchPlain& root); 
  void disorientDebugHelper(TreeAln &traln, const BranchPlain& root); 
  void exa_evaluateGeneric(TreeAln &traln, const BranchPlain& root, bool changedOrientation ); 
  void coreEvalSubTree(TreeAln& traln, const BranchPlain &root); 

private: 			// ATTRIBUTES 
  std::shared_ptr<TreeAln> debugTraln;  
  bool verifyLnl; 
  std::unique_ptr<ArrayPolicy> arrayPolicy; 
  ArrayOrientation arrayOrientation; 
  std::shared_ptr<ArrayReservoir> _arrayReservoir; 
}; 
#endif

