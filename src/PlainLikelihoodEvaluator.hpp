#ifndef _PLAIN_LNL_EVAL
#define _PLAIN_LNL_EVAL

#include "LikelihoodEvaluator.hpp"
#include <cassert>


/** 
    @file PlainLikelihoodEvaluator.hpp

    @brief a likelihood evaluator that does not use any backup arrays  
 */ 
class PlainLikelihoodEvaluator : public LikelihoodEvaluator
{  
public: 
  virtual ~PlainLikelihoodEvaluator(){}

  virtual double evaluate(TreeAln &traln, const BranchPlain &evalBranch,  bool fullTraversal ); 
  virtual double evaluatePartitionsWithRoot( TreeAln &traln, const BranchPlain& root,  const std::vector<nat>& partitions, bool fullTraversal)  {assert(0); return 0; }
  virtual void evalSubtree( TreeAln &traln, const BranchPlain &evalBranch); 

  virtual void imprint(const TreeAln &traln) {} // NOOP 
  virtual void resetToImprinted(TreeAln &traln) {} // NOOP

  virtual void resetSomePartitionsToImprinted(TreeAln &traln, std::vector<nat> partitions) ; 

  virtual std::unique_ptr<LikelihoodEvaluator> clone() const {return std::unique_ptr<LikelihoodEvaluator>(new PlainLikelihoodEvaluator(*this)); }

private: 
  /** 
      @brief does a full traversal and determines which arrays are not
      oriented correctly
      @return indicates whether entire subtree was oriented correctly 
   */
  bool traverseAndTouch(const TreeAln &traln, const BranchPlain &b); 


}; 



#endif
