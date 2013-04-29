#include "CoupledChains.hpp"
#include "Chain.hpp"
#include "randomness.h"
#include "globals.h"
#include "proposals.h"

#include "treeRead.h"

#include "randomTree.h"


CoupledChains::CoupledChains(int seed, int numCoupled, vector<TreeAln*> trees, int runid, initParamStruct *initParams)
  : rand(seed)
{
  assert((nat)numCoupled == trees.size());

  for(int i = 0; numCoupled; ++i)
    {
      Chain *chain = new Chain(rand.generateSeed(),i, runid, trees[i], initParams); 
      chains.push_back(chain);
    }

  // swap info matrix 
  for(int i = 0; i < numCoupled * numCoupled ; ++i)
    swapInfo.push_back(new SuccessCtr());     
}



CoupledChains::~CoupledChains()
{
  for(auto chain : chains) 
    exa_free(chain); 

  for(auto swapI : swapInfo)
    exa_free(swapI); 
}


void CoupledChains::printSwapInfo()
{
  int numCoupledChains = chains.size(); 
  
  int cnt = 0; 
  for(int i = 0; i < numCoupledChains; ++i)
    {
      if(i < numCoupledChains - 1 )
	PRINT("("); 

      for(int j = 0; j < numCoupledChains; ++j)
	{
	  SuccessCtr *ctr = swapInfo[cnt]; 
	  if(i < j )
	    PRINT("%.1f%%,", i,j, 100 * ctr->getRatioOverall());

	  cnt++; 
	}
      if(i < numCoupledChains - 1 )
	PRINT(")"); 
    }
}




void CoupledChains::switchChainState()
{  
  int numChain = chains.size(); 

  if(numChain == 1)
    return;   

  int chainA = drawGlobalRandIntBound(numChain), 
   chainB = chainA; 
  while(chainA == chainB)
    chainB = drawGlobalRandIntBound(numChain); 

  int coupIdA = chains[chainA]->couplingId,
    coupIdB = chains[chainB]->couplingId; 
  

  if(coupIdA > coupIdB)
    swap(coupIdB, coupIdA); 


  /* 
     IMPORTANT TODO

     this currently assumes that we have a non-informative
     prior. Allow for changes!
  */

  double heatA = chains[chainA]->getChainHeat(),
    heatB = chains[chainB]->getChainHeat(); 

  assert(heatA < 1.f || heatB < 1.f); 

  double lnlA = chains[chainA]->traln->getTr()->likelihood,
    lnlB = chains[chainB]->traln->getTr()->likelihood; 

  double 
    aB = lnlA *  heatB,
    bA = lnlB *  heatA,
    aA = lnlA * heatA,
    bB =  lnlB *  heatB; 

  double accRatio = exp(( aB + bA )  - (aA + bB )); 

  Chain *a = chains[ chainA],
    *b = chains[ chainB] ; 
  
  /* do the swap */
  if( drawGlobalDouble01()  < accRatio)
    {
      /* everything, we need to swap */
      swap(a->couplingId, b->couplingId); 

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

      swapInfo[r * chains.size() + c]->accept(); 
            // gAInfo.swapInfo[runId][r * gAInfo.numberCoupledChains + c].accept();
    } 
  else 
    {
      int r = MIN(a->couplingId, b->couplingId ); 
      int c = MAX(a->couplingId, b->couplingId); 

      swapInfo[r * chains.size() + c]->reject(); 

      // gAInfo.swapInfo[runId][r * gAInfo.numberCoupledChains + c].reject();
    }
}


void CoupledChains::executePart(int gensToRun)
{  
  /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
  for(nat i = 0; i < chains.size(); ++i)
    chains[i]->applyChainStateToTree();

  for(int genCtr = 0; genCtr < gensToRun; genCtr += gAInfo.swapInterval)
    {
      boolean timeToTune = FALSE; 

      for(int chainCtr = 0; chainCtr < gAInfo.numberCoupledChains; ++chainCtr)
	{      
	  Chain *curChain = chains[ chainCtr];

	  for(int i = 0; i < gAInfo.swapInterval; ++i)
	    {
	      curChain->step();		  
	      if(gAInfo.tuneFreq > 0 && gAInfo.tuneHeat && curChain->currentGeneration % gAInfo.tuneFreq == gAInfo.tuneFreq - 1 )
		timeToTune = TRUE; 
	    }
	}

      if(timeToTune)
	{	      
	  tuneTemperature();
	}
	    
      
      switchChainState();
    }

  for(int i = 0; i < gAInfo.numberCoupledChains; ++i)
    chains[i]->saveTreeStateToChain();
}


void CoupledChains::tuneTemperature()
{
  /* naive strategy: tune, s.t. the coldest hot chain swaps
     with the coldest chain in 23.4% of all cases */
  SuccessCtr *c = swapInfo[1]; 

  int batch = chains[0]->currentGeneration / gAInfo.tuneFreq; 
  temperature = tuneParameter(batch, c->getRatioInLastInterval(), temperature, FALSE); 

  c->reset();
  
  // update the chains 
  for(auto chain : chains)
    chain->setDeltaT(temperature); 
}




