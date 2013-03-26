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
#include "chain.h"
#include  "adapters.h"
#include "eval.h"
#include "proposals.h"

#include "tune.h"

extern double masterTime; 


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

  double lnlA = chains[chainA].lnl.likelihood,
    lnlB = chains[chainB].lnl.likelihood; 

  double 
    aB = lnlA *  heatB,
    bA = lnlB *  heatA,
    aA = lnlA * heatA,
    bB =  lnlB *  heatB; 

  double accRatio = exp(( aB + bA )  - (aA + bB )); 

  state *a = chains + chainA,
    *b = chains  + chainB ; 
  
  /* do the swap */
  if( drawGlobalDouble01()  < accRatio)
    {

      /* everything, we need to swap */
      swpInt(&(a->couplingId), &(b->couplingId)); 

      for(int i = 0; i < a->numProposals; ++i)
	{
	  proposalFunction *tmp  = a->proposals[i]; 
	  a->proposals[i] = b->proposals[i]; 
	  b->proposals[i] = tmp; 
	}

      int r = MIN(a->couplingId, b->couplingId ); 
      int c = MAX(a->couplingId, b->couplingId); 
      cntAccept(&(gAInfo.swapInfo[runId][r * gAInfo.numberCoupledChains + c ]) );       
    } 
  else 
    {
      int r = MIN(a->couplingId, b->couplingId ); 
      int c = MAX(a->couplingId, b->couplingId); 
      cntReject(&(gAInfo.swapInfo[runId][r * gAInfo.numberCoupledChains + c ]) ); 
    }
}



/**
   @brief Execute a portion of one run. 

   @param chains  -- the pointer to the beginning of the chains in the chain array. 
   @param gensToRun  -- number of generations each chain in this run should proceed. 

 */
void executeOneRun(state *chains, int gensToRun )
{  
  int
    runId  =chains[0].id / gAInfo.numberCoupledChains; 

  if(gAInfo.numberCoupledChains > 1 )
    {
#ifdef MC3_SPACE_FOR_TIME
      /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
      for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
	applyChainStateToTree(chains+i);
#endif
      
      for(int genCtr = 0; genCtr < gensToRun; genCtr += SWITCH_AFTER_GEN)
	{
	  boolean timeToTune = FALSE; 

	  for(int chainCtr = 0; chainCtr < gAInfo.numberCoupledChains; ++chainCtr)
	    {      
	      state *curChain = chains + chainCtr; /* TODO */

#ifndef MC3_SPACE_FOR_TIME
	      applyChainStateToTree(curChain );
#endif
	      for(int i = 0; i < SWITCH_AFTER_GEN; ++i)
		{
		  step(curChain);
		  if(curChain->currentGeneration % TUNE_FREQUENCY == TUNE_FREQUENCY -1 )
		    timeToTune = TRUE; 
		}

#ifndef MC3_SPACE_FOR_TIME
	      saveTreeStateToChain(curChain);
#endif
	    }

	  if(timeToTune)
	    {
	      /* naive strategy: tune, s.t. the coldest hot chain swaps
		 with the coldest chain in 23.4% of all cases */
	      successCtr *c = &(gAInfo.swapInfo[runId][1]); 
	      
	      gAInfo.temperature[runId] = tuneParameter(chains[0].currentGeneration / TUNE_FREQUENCY , getRatioLocal(c), gAInfo.temperature[runId], FALSE);
	      resetCtr(c);
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


