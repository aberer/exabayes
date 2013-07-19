#ifndef _PLAIN_LNL_EVAL
#define _PLAIN_LNL_EVAL

#include "LikelihoodEvaluator.hpp"


/** 
    @file PlainLikelihoodEvaluator.hpp

    @brief a likelihood evaluator that does not use any backup arrays  
 */ 
class PlainLikelihoodEvaluator : public LikelihoodEvaluator
{  
public: 
  virtual double evaluate(TreeAln &traln, const Branch &evalBranch,  bool fullTraversal ); 
  virtual double evaluatePartitions( TreeAln &traln, const std::vector<nat>& partitions); 
  virtual void evalSubtree( TreeAln &traln, const Branch &evalBranch); 

  virtual void imprint(const TreeAln &traln) {} // NOOP 
  virtual void resetToImprinted(TreeAln &traln) {} // NOOP

private: 
  /** 
      @brief does a full traversal and determines which arrays are not
      oriented correctly
      @return indicates whether entire subtree was oriented correctly 
   */
  bool traverseAndTouch(const TreeAln &traln, const Branch &b); 


}; 



#endif
