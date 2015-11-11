#include "ParsimonyEvaluator.hpp"

// #define VERIFY_PARSIMONY

#if HAVE_PLL == 0
#include <mpi.h>
extern MPI_Comm comm; 
#endif

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
#if HAVE_PLL != 0
  newviewParsimony(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p); 
#else 
  newviewParsimony(&(traln.getTrHandle()),p); 
#endif
}


static void verifyParsimony(TreeAln  &traln, nodeptr p, std::array<parsimonyNumber,2> lengths )
{
#ifdef VERIFY_PARSIMONY
  auto partitionParsimony2 = std::array<parsimonyNumber,2>{ { 0,0} };

#if HAVE_PLL == 0 
  evaluateParsimony(&(traln.getTrHandle()),  p , TRUE, partitionParsimony2.data()); 

  MPI_Allreduce(MPI_IN_PLACE, partitionParsimony2.data(), partitionParsimony2.size(), MPI_UNSIGNED, MPI_SUM, comm) ; 
  
#else 
  evaluateParsimony(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p , TRUE, partitionParsimony2.data()); 
#endif

  auto lengthsVerified = partitionParsimony2; 

  if(lengthsVerified[0] != lengths[0] || lengthsVerified[1] != lengths[1])
    {
      std::cout << "error with parsimony: " << lengthsVerified[0] << "," << lengthsVerified[1]<< " != " <<  lengths[0]  << ","  << lengths[1] << std::endl; 
      assert(0); 
    }
#endif
}


auto ParsimonyEvaluator::evaluate(TreeAln &traln, nodeptr p, bool fullTraversal, bool doParallelReduce ) 
  -> std::array<parsimonyNumber,2>
{
  auto stateParsimony = std::array<parsimonyNumber,2>{ {0,0} }; 

#if HAVE_PLL != 0 
  evaluateParsimony(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p,
		    fullTraversal ? TRUE : FALSE , stateParsimony.data());
#else   
  evaluateParsimony(&(traln.getTrHandle()), p, fullTraversal ? TRUE  : FALSE, stateParsimony.data()); 
#endif

  auto result =  stateParsimony; 
#if HAVE_PLL ==  0
  // with examl, we have to reduce 
  if(doParallelReduce)
    {
    MPI_Allreduce(MPI_IN_PLACE, result.data(), result.size(), MPI_UNSIGNED, MPI_SUM, comm);
    // assert(doParallelReduce); 
    }
#endif

  verifyParsimony(traln, p, result); 

  return result; 
}
