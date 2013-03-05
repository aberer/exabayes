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



void mcmc(tree *tr, analdef *adef)
{ 
  
  initRNG(seed);

  /* TODO have removed that -- problematic?   */
  /* assert( isTip(tr->start->number, tr->mxtips )); */

  state *indiChains = NULL; 		/* one state per indipendent run/chain */  
  initParamStruct *initParams = NULL;

  initializeIndependentChains(tr, &indiChains, &initParams); 

  /* TODO  remove this global again */
  numberOfChains = initParams->numIndiChains; 
  int diagFreq = initParams->diagFreq; 

  printf("num indi chains is %d\n", numberOfChains); 

  evaluateGeneric(tr, tr->start, TRUE);
  PRINT( "after reset start: %f\n\n", tr->likelihood );

  boolean hasConverged = FALSE;   
  while(NOT hasConverged)
    {
      for(int i = 0; i < numberOfChains; ++i)
	{
	  state *curChain = indiChains + i;
	  applyChainStateToTree(curChain, tr);

	  for(int j = 0; j < diagFreq; ++j)
	    {
	      step(curChain);
	    }

	  saveTreeStateToChain(curChain, tr);
	}
      
      hasConverged = convergenceDiagnostic(indiChains, numberOfChains); 
    }

  if(processID == 0)
    {
      for(int i = 0; i < numberOfChains; ++i)
	finalizeOutputFiles(indiChains + i);
      PRINT("\nConverged after %d generations\n",  indiChains[0].currentGeneration);
      PRINT("\nTotal execution time: %f seconds\n", gettime() - masterTime); 
    }  

}


