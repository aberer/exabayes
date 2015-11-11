#include "ParsimonyEvaluator.hpp"

void ParsimonyEvaluator::disorientNode(nodeptr p)
{
  if(p->xPars == 1 )
    {
      p->xPars = 0; 
      p->next->xPars = 1 ; 
      p->next->next->xPars = 0 ; 
    }
}


void ParsimonyEvaluator::evaluateSubtree(TreeAln &traln, nodeptr p)
{
  newviewParsimony(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p); 
}


std::array<nat,2> ParsimonyEvaluator::evaluate(TreeAln &traln, nodeptr p, bool fullTraversal ) 
{
  auto result = std::array<parsimonyNumber,2>{{0,0}}; 

  auto &pHandle = traln.getPartitionsHandle();
  auto &tHandle = traln.getTrHandle();
  
  assert(not tHandle.fastParsimony ); 

  evaluateParsimony(&tHandle, &pHandle, p, fullTraversal ? PLL_TRUE : PLL_FALSE);

  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      auto &partition = traln.getPartition(i); 
      // HACK 
      switch(partition.getStates())
	{
	case 4: 
	  result[0] += partition.getPartitionParsimony(); 
	  break; 
	case 20: 
	  result[1] += partition.getPartitionParsimony(); 
	  break; 
	default : 
	  assert(0); 
	}
    }

  // meh 
  // if(doParallelReduce)
  //   {
  //     auto tmp = std::vector<nat> {result[0], result[1]}; 
  //     _commPtr->allReduce(tmp); 
  //     result[0]= tmp[0]; 
  //     result[1] = tmp[1]; 
  //   }
  

  return result; 
}
