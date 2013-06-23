#include "ParsimonyEvaluator.hpp"

void ParsimonyEvaluator::evaluateSubtree(TreeAln &traln, nodeptr p)
{
#if HAVE_PLL != 0
  newviewParsimony(traln.getTr(), traln.getPartitionsPtr(), p); 
#else 
  newviewParsimony(traln.getTr(),p); 
#endif
}
 
void ParsimonyEvaluator::evaluate(TreeAln &traln, nodeptr p, bool fullTraversal, vector<nat> &partitionParsimony)
{
  partitionParsimony.clear();   
  nat* tmp = (nat*)exa_calloc(traln.getNumberOfPartitions(), sizeof(nat)); 
  
#if HAVE_PLL != 0 
  evaluateParsimony(traln.getTr(), traln.getPartitionsPtr(), p, fullTraversal ? TRUE : FALSE , tmp); 
#else 
  evaluateParsimony(traln.getTr(), p, fullTraversal ? TRUE  : FALSE, tmp); 
  // assert(0);
#endif

  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    partitionParsimony.push_back(tmp[i]); 

  exa_free(tmp); 
}
