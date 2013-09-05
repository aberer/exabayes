#ifndef PARSIMONY_EVALUATOR
#define PARSIMONY_EVALUATOR

#include "TreeAln.hpp"

class ParsimonyEvaluator
{
public:
  void evaluateSubtree(TreeAln &traln, nodeptr p); 

  /** 
      @brief evaluates the parsimony score of the tree

      @param parsimonyLength the per partition parsimony score for the
      transition between the descendent nodes
   */ 
  void evaluate(TreeAln &traln, nodeptr p, bool fullTraversal, std::vector<nat> &partitionParsimony, std::vector<nat> &parsimonyLength); 

private: 


}; 

#endif
