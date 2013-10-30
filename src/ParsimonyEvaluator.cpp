#include "ParsimonyEvaluator.hpp"

void ParsimonyEvaluator::evaluateSubtree(TreeAln &traln, nodeptr p)
{
#if HAVE_PLL != 0
  newviewParsimony(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p); 
#else 
  newviewParsimony(&(traln.getTrHandle()),p); 
#endif
}
 
void ParsimonyEvaluator::evaluate(TreeAln &traln, nodeptr p, bool fullTraversal, 
				  std::vector<nat> &partitionParsimony, std::vector<nat> &parsimonyLength)
{
  partitionParsimony = std::vector<nat>(traln.getNumberOfPartitions(), 0);
  parsimonyLength = std::vector<nat>(traln.getNumberOfPartitions(),0);

#if HAVE_PLL != 0 
  evaluateParsimony(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p,
		    fullTraversal ? TRUE : FALSE , 
		    &(partitionParsimony[0]), &(parsimonyLength[0]));
#else   
  evaluateParsimony(&(traln.getTrHandle()), p, fullTraversal ? TRUE  : FALSE, &(partitionParsimony[0]), &(parsimonyLength[0])); 
#endif

}
