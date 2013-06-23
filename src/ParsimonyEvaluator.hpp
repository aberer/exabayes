#ifndef PARSIMONY_EVALUATOR

#include "TreeAln.hpp"

class ParsimonyEvaluator
{
public:
  void evaluateSubtree(TreeAln &traln, nodeptr p); 
  void evaluate(TreeAln &traln, nodeptr p, bool fullTraversal, vector<nat> &partitionParsimony);   

private: 


}; 

#endif