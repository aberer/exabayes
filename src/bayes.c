#include "axml.h"
#include "proposalStructs.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"

#include "proposals.h"
#include "convergence.h"

#include "bayes-topo.h"


/* TODO outsource  */
#include "chain.h"

#include  "adapterCode.h"

#include "eval.h"

extern double masterTime; 


#ifdef _USE_NCL_PARSER
#include "nclConfigReader.h"
void addInitParameters(state *curstate, initParamStruct *initParams)
{
  curstate->proposalWeights[E_SPR] = initParams->initSPRWeight; 
  curstate->proposalWeights[UPDATE_GAMMA] = initParams->initGammaWeight;
  curstate->proposalWeights[UPDATE_GAMMA_EXP] = initParams->initGammaExpWeight;
  curstate->proposalWeights[UPDATE_MODEL] = initParams->initModelWeight; 
  curstate->proposalWeights[UPDATE_SINGLE_BL] = initParams->initSingleBranchWeight; 
  curstate->proposalWeights[UPDATE_SINGLE_BL_EXP] = initParams->initSingleBranchExpWeight; 
curstate->proposalWeights[UPDATE_SINGLE_BL_BIUNIF] = initParams->initSingleBranchBiunifWeight;
curstate->proposalWeights[UPDATE_MODEL_BIUNIF] = initParams->initModelBiunifWeight;
curstate->proposalWeights[UPDATE_MODEL_SINGLE_BIUNIF] = initParams->initModelSingleBiunifWeight;
curstate->proposalWeights[UPDATE_MODEL_ALL_BIUNIF] = initParams->initModelAllBiunifWeight;
curstate->proposalWeights[UPDATE_MODEL_PERM_BIUNIF] = initParams->initModelPermBiunifWeight;
curstate->proposalWeights[UPDATE_FREQUENCIES_BIUNIF] = initParams->initFrequenciesWeight;
curstate->proposalWeights[E_SPR_MAPPED] = initParams->initEsprMappedWeight;
  //PROPOSALADD addInitParameters NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
  curstate->numGen = initParams->numGen; 
  curstate->penaltyFactor = initParams->initPenaltyFactor; 
  curstate->samplingFrequency = initParams->samplingFrequency; 
  curstate->eSprStopProb = initParams->eSprStopProb; 
}
#else 
int parseConfig(state *theState);
#endif



void initDefaultValues(state *theState, tree *tr)
{
  theState->curprior = 1; 
  theState->hastings = 1; 
  theState->currentGeneration = 0; 

  theState->brLenRemem.bl_sliding_window_w = 0.005;
  theState->brLenRemem.bl_prior = 1.0;
  theState->brLenRemem.bl_prior_exp_lambda = 0.1 ;
  //this can be extended to more than one partition, but one for now
  
  theState->modelRemem.model = 0;
  theState->modelRemem.rt_sliding_window_w = 0.5;

  pInfo *partition = getPartition(theState,theState->modelRemem.model); 

  theState->modelRemem.nstates = partition->states; /* 4 for DNA */
  theState->modelRemem.numSubsRates = (theState->modelRemem.nstates * theState->modelRemem.nstates - theState->modelRemem.nstates) / 2; /* 6 for DNA */
  theState->modelRemem.curSubsRates = (double *) exa_malloc(theState->modelRemem.numSubsRates * sizeof(double));
  
  theState->frequRemem.model = 0;
  theState->frequRemem.numFrequRates = partition->states; /* 4 for DNA */
  theState->frequRemem.curFrequRates = (double *) exa_malloc(theState->frequRemem.numFrequRates * sizeof(double));
  
  theState->gammaRemem.gm_sliding_window_w = 0.75;
  theState->brLenRemem.single_bl_branch = -1;


  theState->proposalWeights[E_SPR] = 0.0; 
  theState->proposalWeights[UPDATE_MODEL] = 0.0; 
  theState->proposalWeights[UPDATE_GAMMA] = 0.0; 
  theState->proposalWeights[UPDATE_GAMMA_EXP] = 0.0; 
  theState->proposalWeights[UPDATE_SINGLE_BL] = 0.0;   
  theState->proposalWeights[UPDATE_SINGLE_BL_EXP] = 0.0;   
theState->proposalWeights[UPDATE_SINGLE_BL_BIUNIF] = 0.0;
theState->proposalWeights[UPDATE_MODEL_BIUNIF] = 0.0;
theState->proposalWeights[UPDATE_MODEL_SINGLE_BIUNIF] = 0.0;
theState->proposalWeights[UPDATE_MODEL_ALL_BIUNIF] = 0.0;
theState->proposalWeights[UPDATE_MODEL_PERM_BIUNIF] = 0.0;
theState->proposalWeights[UPDATE_FREQUENCIES_BIUNIF] = 0.0;
theState->proposalWeights[E_SPR_MAPPED] = 0.0;
  //PROPOSALADD initDefaultValues NOTE Do not remove/modify  this line. The script addProposal.pl needs it as an identifier.
  
  theState->numGen = 1000000;
  theState->penaltyFactor = 0.0;
}




/* TODO adapt likelihood */
void switchChainState(state *chains)
{
  int numChain = gAInfo.numberCoupledChains; 

  if(numChain == 1)
    return;   
  
  /* randCtr_t r = drawGlobalRandInt();  */
    
  int chainA = drawGlobalRandIntBound(numChain), 
   chainB = chainA; 
  while(chainA == chainB)
    chainB = drawGlobalRandIntBound(numChain); 

  /* int absIdA =  chains[chainA].id,  */
  /*   absIdB = chains[chainB].id;  */

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

  double accRatio = ( aB + bA )  - (aA + bB ); 

  /* do the swap */
  if( drawGlobalDouble01()  < exp(accRatio))
    {

      int tmp = chains[chainA].couplingId ;
      chains[chainA].couplingId =  chains[chainB].couplingId;
      chains[chainB].couplingId = tmp ;

      int tmpArray[NUM_PROPOSALS] ; 
      memcpy(tmpArray, chains[chainB].acceptedProposals, sizeof(int) * NUM_PROPOSALS); 
      memcpy(chains[chainB].acceptedProposals, chains[chainA].acceptedProposals, sizeof(int) * NUM_PROPOSALS); 
      memcpy(chains[chainA].acceptedProposals, tmpArray, sizeof(int) * NUM_PROPOSALS); 

      memcpy(tmpArray, chains[chainB].rejectedProposals, sizeof(int) * NUM_PROPOSALS); 
      memcpy(chains[chainB].rejectedProposals, chains[chainA].rejectedProposals, sizeof(int) * NUM_PROPOSALS); 
      memcpy(chains[chainA].rejectedProposals, tmpArray, sizeof(int) * NUM_PROPOSALS); 

      printf("%d switches with %d\n", chainA, chainB); 

/* #if HAVE_PLL == 1        */
/*       partitionList *pTmp = gAInfo.partitions[chainA];  */
/*       gAInfo.partitions[chainA] = gAInfo.partitions[chainB]; */
/*       gAInfo.partitions[chainB] = pTmp;       */
/* #endif */

      gAInfo.successFullSwitchesBatch++; 
    } 
}




/**
   @brief Execute a portion of one run. 

   @param state *chains  -- the pointer to the beginning of the chains in the chain array. 
   @param int gensToRun  -- number of generations each chain in this run should proceed. 

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



void runChains(state *allChains,   int diagFreq)
{
  boolean hasConverged = FALSE;   
  while(NOT hasConverged)
    {      
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	executeOneRun(allChains + (i * gAInfo.numberCoupledChains), diagFreq); 
      
      hasConverged = convergenceDiagnostic(allChains, gAInfo.numberOfRuns); 
    }
}


void exa_main(tree *tr, analdef *adef)
{   
  state *indiChains = NULL; 		/* one state per indipendent run/chain */  
  initParamStruct *initParams = NULL;

  timeIncrement = gettime();
  
  initializeIndependentChains(tr, adef,  &indiChains, &initParams); 

  assert(gAInfo.numberCoupledChains > 0);
  /* TODO more bla bla  */  

  int diagFreq = initParams->diagFreq; 
  
  gAInfo.allChains = indiChains; 

  runChains(indiChains, diagFreq); 

  if(processID == 0)
    {
      for(int i = 0; i < gAInfo.numberOfRuns; ++i)
	finalizeOutputFiles(indiChains + i);
      PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    } 
}


