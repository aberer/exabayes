#ifndef PARSIMONY_EVALUATOR
#define PARSIMONY_EVALUATOR

#include "comm/RemoteComm.hpp"
#include <cassert>
#include "model/TreeAln.hpp"

class ParsimonyEvaluator
{
public:

  void evaluateSubtree(TreeAln &traln, nodeptr p);

  /** 
      @brief evaluates the parsimony score of the tree

      @param parsimonyLength the per partition parsimony score for the
      transition between the descendent nodes
   */ 
  auto evaluate(TreeAln &traln, nodeptr p, bool fullTraversal )   
    -> std::array<parsimonyNumber,2>;

  static nat numState2pos(nat numState) 
  { 
    switch(numState)
      {
      case 4: 
	return 0 ; 
      case 20: 
	return 1; 
      default: 
	assert(0); 
      }
  }

  static  nat pos2numstate(nat pos)
  {
    switch(pos)
      {
      case 0: 
	return 4; 
      case 1: 
	return 20; 
      default : 
	assert(0); 
      }
  }

  static void disorientNode(nodeptr p); 

private: 
  std::shared_ptr<RemoteComm> _commPtr; 
}; 

#endif
