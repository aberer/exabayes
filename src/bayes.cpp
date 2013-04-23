/**
   @file bayes.c
 
   @brief The top-level functionality of ExaBayes.

*/


#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"
#include "proposals.h"
#include "nclConfigReader.h"
#include "misc-utils.h"
#include "chain.h"
#include "adapters.h"
#include "eval.h"
#include "proposals.h"
#include "tune.h"
#include "prsfComputer.h"
#include "TreeAln.hpp"
#include "AvgSplitFreqAssessor.hpp"

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

  double lnlA = chains[chainA].traln->getTr()->likelihood,
    lnlB = chains[chainB].traln->getTr()->likelihood; 

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

      FILE *tmp = a-> topologyFile; 
      a->topologyFile = b->topologyFile; 
      b->topologyFile = tmp; 
      tmp = a->outputParamFile; 
      a->outputParamFile = b->outputParamFile; 
      b->outputParamFile = tmp; 

      int r = MIN(a->couplingId, b->couplingId ); 
      int c = MAX(a->couplingId, b->couplingId); 
      
      gAInfo.swapInfo[runId][r * gAInfo.numberCoupledChains + c].accept();
    } 
  else 
    {
      int r = MIN(a->couplingId, b->couplingId ); 
      int c = MAX(a->couplingId, b->couplingId); 

      gAInfo.swapInfo[runId][r * gAInfo.numberCoupledChains + c].reject();
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
      /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
      for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
	applyChainStateToTree(chains+i);
      
      for(int genCtr = 0; genCtr < gensToRun; genCtr += gAInfo.swapInterval)
	{
	  boolean timeToTune = FALSE; 

	  for(int chainCtr = 0; chainCtr < gAInfo.numberCoupledChains; ++chainCtr)
	    {      
	      state *curChain = chains + chainCtr; /* TODO */

	      for(int i = 0; i < gAInfo.swapInterval; ++i)
		{
		  step(curChain);		  
		  if(gAInfo.tuneFreq > 0 && gAInfo.tuneHeat && curChain->currentGeneration % gAInfo.tuneFreq == gAInfo.tuneFreq - 1 )
		    timeToTune = TRUE; 
		}
	    }

	  if(timeToTune)
	    {	      
	      /* naive strategy: tune, s.t. the coldest hot chain swaps
		 with the coldest chain in 23.4% of all cases */
	      SuccessCtr *c = &(gAInfo.swapInfo[runId][1]); 
	      
	      gAInfo.temperature[runId] = tuneParameter(chains[0].currentGeneration / gAInfo.tuneFreq , c->getRatioInLastInterval(), gAInfo.temperature[runId], FALSE);
	      c->reset();
	    }
	    
	  switchChainState(chains);
	}

      for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
	saveTreeStateToChain(chains+i);

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





#include <sstream>
bool convergenceDiagnostic(state *allChains)
{
  if(gAInfo.numberOfRuns > 1)
    { 
      vector<string> fns; 
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	{
	  stringstream ss; 
	  ss <<  PROGRAM_NAME << "_topologies." << run_id << "." << i; 
	  fns.push_back(ss.str());
	}
     
      AvgSplitFreqAssessor asdsf(fns);

      int end = asdsf.getEnd();
      
      int treesInBatch = gAInfo.diagFreq / gAInfo.samplingFrequency; 

      end /= treesInBatch; 
      end *= treesInBatch;       

      if(end > 0)
	{	  
	  asdsf.setEnd(end);
	  if( gAInfo.burninGen > 0 )
	    {
	      assert(gAInfo.burninProportion == 0.); 

	      int treesToDiscard =  gAInfo.burninGen / gAInfo.samplingFrequency; 

	      if(end < treesToDiscard + 2 )
		return false; 
	      else 
		asdsf.setStart(treesToDiscard);  
	    }
	  else 
	    {
	      assert(gAInfo.burninGen == 0); 
	      int start = (int)((double)end * gAInfo.burninProportion  ); 
	      asdsf.setStart(start);
	    } 

	  asdsf.extractBips();
	  double asdsfVal = asdsf.computeAsdsf(gAInfo.asdsfIgnoreFreq);
      
	  if(processID == 0)
	    cout << "ASDSF for trees " << asdsf.getStart() << "-" << asdsf.getEnd() << ": " << asdsfVal << endl; 

	  return asdsfVal < gAInfo.asdsfConvergence; 

	}
      else 
	return false; 
      
    }
  else 
    return allChains[0].currentGeneration > gAInfo.numGen; 
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

      hasConverged = convergenceDiagnostic(allChains); 

#ifdef ENABLE_PRSF
      if(isOutputProcess)
	printPRSF(run_id);
#endif
    }
}



// #define TEST 
 
#include "MyTestProposal.hpp"



// /**
//    @brief the main ExaBayes function.

//    @param tr -- a tree structure that has been initialize in one of the adapter mains. 
//    @param adef -- the legacy adef
//  */
void exa_main (analdef *adef, int seed)
{   
  state *indiChains = NULL; 		/* one state per indipendent run/chain */  

  timeIncrement = gettime();
  gAInfo.adef = adef; 

#ifdef TEST   
  MyTestProposal prop();
  
  

  exit(0); 
#endif


  initializeIndependentChains(adef,  seed, &indiChains); 

  assert(gAInfo.numberCoupledChains > 0);

  gAInfo.allChains = indiChains; 

  assert(gAInfo.diagFreq != 0 ); 
  runChains(indiChains, gAInfo.diagFreq); 

  if(isOutputProcess() )
    {
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	finalizeOutputFiles(indiChains + i);
      PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    } 
}


