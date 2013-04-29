#include <sstream>


#include "CoupledChains.hpp"
#include "Chain.hpp"
#include "randomness.h"
#include "globals.h"
#include "proposals.h"
#include "treeRead.h"
#include "output.h"
#include "topology-utils.h"



CoupledChains::CoupledChains(int seed, int numCoupled, vector<TreeAln*> trees, int runid, initParamStruct *initParams)
  : rand(seed)
{
  assert((nat)numCoupled == trees.size());

  for(int i = 0; i < numCoupled; ++i)
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

  int chainA = rand.drawRandInt(numChain),
    chainB = chainA; 
  while(chainA == chainB)
    chainB = rand.drawRandInt(numChain); 

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

  assert(heatA <= 1.f || heatB <= 1.f); 

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
  if( rand.drawRandDouble01()  < accRatio)
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


void CoupledChains::chainInfo()
{
  // find cold chain
  Chain *coldChain = NULL; 
  for(auto chain : chains )
    if(chain->couplingId == 0)
      coldChain = chain; 
  assert(coldChain != NULL); 
  
  tree *tr = coldChain->traln->getTr(); 

  PRINT( "[run: %d] [time %.2f] gen: %d Likelihood: %.2f\tTL=%.2f\t",runid,   gettime()  - timeIncrement  , coldChain->currentGeneration, coldChain->traln->getTr()->likelihood, branchLengthToReal(tr, getTreeLength(coldChain->traln, tr->nodep[1]->back)));

  
  // print hot chains
  vector<Chain*> sortedChains(chains.size()); 
  for(auto chain : chains)
    sortedChains[chain->couplingId] = chain; 
  for(nat i = 1 ; i < chains.size(); ++i)
    {
      Chain *chain = sortedChains[i]; 
      double heat = chain->getChainHeat();
      assert(heat < 1.0f); 
      
      PRINT("lnl_beta(%.2f)=%.2f\t", heat, chain->traln->getTr()->likelihood); 
    }

  PRINT("\n"); 

  /* just output how much time has passed since the last increment */
  timeIncrement = gettime(); 	

  char* names[]= {"TOPO", "BL", "FREQ", "MODEL", "HETER"} ;   

  for(int i = 1; i < NUM_PROP_CATS + 1 ; ++i)
    {
      category_t type = category_t(i); 
      PRINT("%s:", names[i-1]);       
      for(int j = 0; j < coldChain->numProposals;++j )
	{
	  proposalFunction *pf = coldChain->proposals[j]; 	
	  if(type == pf->category)
	    {
	      stringstream buf; 
	      buf << pf->sCtr; 
	      PRINT("\t%s: %s\t", pf->name, buf.str().c_str()); 
	    }
	}
      PRINT("\n"); 
    }

  PRINT("\n");
}


void CoupledChains::executePart(int gensToRun)
{  
  /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
  for(nat i = 0; i < chains.size(); ++i)
    chains[i]->applyChainStateToTree();

  for(int genCtr = 0; genCtr < gensToRun; genCtr += gAInfo.swapInterval)
    {
      bool timeToTune = false; 
      bool timeToPrint = false; 

      for(int chainCtr = 0; chainCtr < gAInfo.numberCoupledChains; ++chainCtr)
	{      
	  Chain *chain = getChain(chainCtr); 

	  for(int i = 0; i < gAInfo.swapInterval; ++i)
	    {
	      chain->step(); 
	      timeToTune = gAInfo.tuneFreq > 0 && gAInfo.tuneHeat && chain->currentGeneration % gAInfo.tuneFreq == gAInfo.tuneFreq - 1 ; 
	      timeToPrint = isOutputProcess() && chain->couplingId == 0 && gAInfo.printFreq > 0 && chain->currentGeneration % gAInfo.printFreq == gAInfo.printFreq - 1 ; // that is a bit nonsensical 
	    }
	}

      if(timeToTune)
	tuneTemperature();
      
      if(timeToPrint)
	
	    
      
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




