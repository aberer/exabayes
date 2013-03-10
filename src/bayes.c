#include "common.h"

#include "axml.h"
#include "proposalStructs.h"
#include "randomness.h"
#include "globals.h"
#include "main-common.h"
#include "output.h"

#include "proposals.h"
#include "convergence.h"


/* TODO outsource  */
#include "chain.h"


extern double masterTime; 


void saveTreeStateToChain(state *chain, tree *tr); 
void applyChainStateToTree(state *chain, tree *tr); 

/* TODO commented this out, since there are some problems with it,
   when we build with the PLL */
/* #define WITH_PERFORMANCE_MEASUREMENTS */



/* call this for verification after the lnl has been evaluated somehow */
void expensiveVerify(tree *tr)
{
#ifdef DEBUG_LNL_VERIFY
  double val1 = tr->likelihood; 
  evaluateGeneric(tr, tr->start, TRUE); 

  if(processID == 0)
    {
      if(fabs (tr->likelihood - val1 ) > 0.1)
      printf("WARNING: found in expensive evaluation: likelihood difference is %f (with before/after)\t%f\t%f\n", fabs (tr->likelihood - val1 ), val1, tr->likelihood); 
      assert(fabs (tr->likelihood - val1 ) < 0.1);   
    }
  
#endif
}




void makeRandomTree(tree *tr); 

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
  theState->modelRemem.nstates = tr->partitionData[theState->modelRemem.model].states; /* 4 for DNA */
  theState->modelRemem.numSubsRates = (theState->modelRemem.nstates * theState->modelRemem.nstates - theState->modelRemem.nstates) / 2; /* 6 for DNA */
  theState->modelRemem.curSubsRates = (double *) malloc(theState->modelRemem.numSubsRates * sizeof(double));
  
  theState->frequRemem.model = 0;
  theState->frequRemem.numFrequRates =tr->partitionData[theState->frequRemem.model].states; /* 4 for DNA */
  theState->frequRemem.curFrequRates = (double *) malloc(theState->frequRemem.numFrequRates * sizeof(double));
  
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
void switchChainState(state *chains, int numChain)
{
  if(numChain == 1)
    return;   

  int chainA = drawRandInt(numChain); 
  int chainB = chainA; 
  while(chainA == chainB)
    chainB = drawRandInt(numChain); 

  /* 
     IMPORTANT TODO

     this currently assumes that we have a non-informative
     prior. Allow for changes!
  */

  double heatA = getChainHeat(chains + chainA ) , 
    heatB  = getChainHeat(chains + chainB); 
  
  double lnlA = chains[chainA].likelihood,
    lnlB = chains[chainB].likelihood; 
  
  double 
    aB = lnlA *  heatB,
    bA = lnlB *  heatA,
    aA = lnlA * heatA,
    bB =  lnlB *  heatB; 

  double accRatio = ( aB + bA )  - (aA + bB ); 

  /* printf("%f,%f,%f,%f\t%f,%f,%f,%f\taccRatio = %f\n", lnlA, lnlB, heatA, heatB,  aB, bA, aA, bB, accRatio);  */

  /* do the swap */
  if( drawRandDouble()  < accRatio)
    {
      int tmp = chains[chainA].couplingId ; 
      chains[chainA].couplingId =  chains[chainB].couplingId; 
      chains[chainB].couplingId = tmp ; 
      
      /* if(processID == 0) */
      /* 	printf("coupled chains %d  and %d switch\n", chainA, chainB);  */
    } 
}



void executeOneRun(state *chains, int gensToRun )
{
  tree *tr = chains[0].tr; 

  for(int genCtr = 0; genCtr < gensToRun; genCtr += SWITCH_AFTER_GEN)
    {
      for(int chainCtr = 0; chainCtr < numberCoupledChains; ++chainCtr)
	{      
	  state *curChain = chains + chainCtr; /* TODO  */
	  applyChainStateToTree(curChain, tr);

	  for(int i = 0; i < SWITCH_AFTER_GEN; ++i)
	    step(curChain);
	  	  
	  saveTreeStateToChain(curChain, tr);
	}

      switchChainState(chains, numberCoupledChains);
    }
}




void mcmc(tree *tr, analdef *adef)
{ 
  
  initRNG(seed);

  /* TODO have removed that -- problematic?   */
  /* assert( isTip(tr->start->number, tr->mxtips )); */

  state *indiChains = NULL; 		/* one state per indipendent run/chain */  
  initParamStruct *initParams = NULL;

  
  initializeIndependentChains(tr, &indiChains, &initParams); 

  /* TODO  remove this global again */
  numberOfRuns = initParams->numIndiChains; 
  numberCoupledChains = initParams->numCoupledChains; 

  int diagFreq = initParams->diagFreq; 

  boolean hasConverged = FALSE;   
  while(NOT hasConverged)
    {      
      for(int i = 0; i < numberOfRuns; ++i)
	executeOneRun(indiChains + (i * numberCoupledChains), diagFreq); 
      
      hasConverged = convergenceDiagnostic(indiChains, numberOfRuns); 
    }

  if(processID == 0)
    {
      for(int i = 0; i < numberOfRuns; ++i)
	finalizeOutputFiles(indiChains + i);
      PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    } 
}


