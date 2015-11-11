#include "ParsimonyEvaluator.hpp"

void ParsimonyEvaluator::evaluateSubtree(TreeAln &traln, nodeptr p)
{
#if HAVE_PLL != 0
  newviewParsimony(traln.getTr(), traln.getPartitionsPtr(), p); 
#else 
  newviewParsimony(traln.getTr(),p); 
#endif
}
 
void ParsimonyEvaluator::evaluate(TreeAln &traln, nodeptr p, bool fullTraversal, std::vector<nat> &partitionParsimony)
{
  partitionParsimony.clear();   
  partitionParsimony.insert(partitionParsimony.begin(), traln.getNumberOfPartitions(), 0); 
  
#if HAVE_PLL != 0 
  evaluateParsimony(traln.getTr(), traln.getPartitionsPtr(), p, fullTraversal ? TRUE : FALSE , &(partitionParsimony[0])); 
#else 
  evaluateParsimony(traln.getTr(), p, fullTraversal ? TRUE  : FALSE, &(partitionParsimony[0])); 
#endif

}
