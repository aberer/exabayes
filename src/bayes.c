/**
   @file bayes.c
 
   @brief The top-level  functionality of ExaBayes.

*/


#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"
#include "proposals.h"
#include "convergence.h"
#include "exa-topology.h"	
#include "nclConfigReader.h"
#include "misc-utils.h"

/* TODO outsource  */
#include "chain.h"

#include  "adapters.h"
#include "eval.h"

#include "proposals.h"

extern double masterTime; 





/* int mapToTriangularIndex(int row, int col) */
/* { */
/*   for(int i = 0; i < row ; ++i) */
/*     for(int j = 0; hp) */
/* } */


/**
   @brief attempt a MC3 switch of chains 

   @param chains -- the pointer to the first chain struct in the chain-array.   
   
   NOTICE does not work with priors yet 
 */
void switchChainState(state *chains)
{  
  int runId = chains[0].id / gAInfo.numberCoupledChains; 

  int numChain = gAInfo.numberCoupledChains; 

  if(numChain == 1)
    return;   

  int chainA = drawGlobalRandIntBound(numChain), 
   chainB = chainA; 
  while(chainA == chainB)
    chainB = drawGlobalRandIntBound(numChain); 

  int coupIdA = chains[chainA].couplingId,
    coupIdB = chains[chainB].couplingId; 

  if(coupIdA > coupIdB)
    swpInt(&coupIdB, &coupIdA); 

  /* 
     IMPORTANT TODO

     this currently assumes that we have a non-informative
     prior. Allow for changes!
  */
  
  double heatA = getChainHeat(chains + chainA ) , 
    heatB  = getChainHeat(chains + chainB); /*  */

  assert(heatA < 1.f || heatB < 1.f); 

  double lnlA = chains[chainA].likelihood,
    lnlB = chains[chainB].likelihood; 

  double 
    aB = lnlA *  heatB,
    bA = lnlB *  heatA,
    aA = lnlA * heatA,
    bB =  lnlB *  heatB; 

  double accRatio = exp(( aB + bA )  - (aA + bB )); 

  /* do the swap */
  if( drawGlobalDouble01()  < accRatio)
    {
      state *a = chains + chainA,
	*b = chains  + chainB ; 

      swpInt(&(a->couplingId), &(b->couplingId)); 
      for(int i = 0; i < a->numProposals;++i)
	{
	  swpInt(&(a->proposals[i]->successCtr.acc), &(b->proposals[i]->successCtr.acc)); 
	  swpInt(&(a->proposals[i]->successCtr.rej), &(b->proposals[i]->successCtr.rej)); 
	}

/*      =0 1 2 3  */
/* ================ */
/*    0 =  0 1 2 */
/*    1 =    3 4  */
/*    2 =      5  */

      /* 	int n = gAInfo.numberCoupledChains;    */
      /* return row* n - (row-1)* row/2 + col - row;   */

      /* TODO */
      /* gAInfo.swapInfo[] */
    } 
  else 
    {
      /* TODO  */
    }
}



/**
   @brief Execute a portion of one run. 

   @param chains  -- the pointer to the beginning of the chains in the chain array. 
   @param gensToRun  -- number of generations each chain in this run should proceed. 

 */
void executeOneRun(state *chains, int gensToRun )
{
  if(gAInfo.numberCoupledChains > 1 )
    {
#ifdef MC3_SPACE_FOR_TIME
      /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
      for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
	applyChainStateToTree(chains+i);
#endif


      for(int genCtr = 0; genCtr < gensToRun; genCtr += SWITCH_AFTER_GEN)
	{
	  for(int chainCtr = 0; chainCtr < gAInfo.numberCoupledChains; ++chainCtr)
	    {      
	      state *curChain = chains + chainCtr; /* TODO */

#ifndef MC3_SPACE_FOR_TIME
	      applyChainStateToTree(curChain );
#endif

	      for(int i = 0; i < SWITCH_AFTER_GEN; ++i)
		step(curChain);
	  	  
#ifndef MC3_SPACE_FOR_TIME
	      saveTreeStateToChain(curChain);
#endif
	    }

	  switchChainState(chains);
	}


#ifdef MC3_SPACE_FOR_TIME
      for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
	saveTreeStateToChain(chains+i);
#endif

    }
  else 
    {      
      state *curChain = chains; 
      applyChainStateToTree(curChain );
      
      for(int genCtr = 0; genCtr < gensToRun; genCtr++)
	step(curChain);

      saveTreeStateToChain(curChain); 
    }
}



/**
   @brief run all chains for a given number of generations
 */
void runChains(state *allChains, int diagFreq)
{
  boolean hasConverged = FALSE;   
  while(NOT hasConverged)
    {      
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	{
	  state *relChains =  allChains + (i * gAInfo.numberCoupledChains); 
	  executeOneRun(relChains, diagFreq); 
	  
	  for(int j = 0; j < gAInfo.numberCoupledChains; ++j)
	    resetSuccessCounters(relChains + j); 
	}

      hasConverged = convergenceDiagnostic(allChains, gAInfo.numberOfRuns); 
    }
}



/**
   @brief the main ExaBayes function.

   @param tr -- a tree structure that has been initialize in one of the adapter mains. 
   @param adef -- the legacy adef
 */
void exa_main(tree *tr, analdef *adef)
{   
  state *indiChains = NULL; 		/* one state per indipendent run/chain */  

  timeIncrement = gettime();
  
  initializeIndependentChains(tr, adef,  &indiChains); 

  assert(gAInfo.numberCoupledChains > 0);
  /* TODO more bla bla  */  

  gAInfo.allChains = indiChains; 
  
  assert(gAInfo.diagFreq != 0 ); 
  runChains(indiChains, gAInfo.diagFreq); 

  if(processID == 0)
    {
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	finalizeOutputFiles(indiChains + i);
      PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    } 
}


