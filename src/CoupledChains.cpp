#include <sstream>

#include "CoupledChains.hpp"
#include "Chain.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "treeRead.h"
#include "AbstractProposal.hpp"
#include "PriorBelief.hpp"

#include "time.hpp"

CoupledChains::CoupledChains(randCtr_t seed, int runNum, string workingdir, int numCoupled,  vector<Chain> &_chains )
  : chains(_chains)
  , heatIncrement(0.1) 
  , rand(seed)
  , runid(runNum) 
  , tuneHeat(false)
  , printFreq(500)
  , swapInterval(1)
  , samplingFreq(100)
  , runname("standardId") 
  , numCoupled(numCoupled) 
{
  // swap info matrix 
  for(int i = 0; i < numCoupled * numCoupled ; ++i)
    swapInfo.push_back(new SuccessCounter());       
  
  stringstream tNameBuilder,
    pNameBuilder;
  tNameBuilder << workingdir << PROGRAM_NAME << "_topologies." << runname << "."  << runid ; 
  pNameBuilder << workingdir << PROGRAM_NAME << "_parameters." << runname << "."  << runid  ; 
  
  /* todo binary Chain file?  */
  topoFile = fopen(tNameBuilder.str().c_str(), "w"); 
  paramFile = fopen(pNameBuilder.str().c_str(), "w");   

  // TODO 

  chains[0].printNexusTreeFileStart(topoFile);
  chains[0].printParamFileStart(paramFile) ;

  // for(auto &c : getChains())
  //   {
  //     cout << c.getTraln() << endl; 
  //   }

}

void CoupledChains::seedChains()
{
  for(auto &c : chains)
    c.reseed(rand.generateSeed()); 
}


CoupledChains::~CoupledChains()
{
}


void CoupledChains::printSwapInfo()
{
  if(chains.size( )== 1 )
    return ; 
  
  int numCoupledChains = chains.size(); 
  
  int cnt = 0; 
  for(int i = 0; i < numCoupledChains; ++i)
    {
      bool isFirst = true; 

      if(i < numCoupledChains - 1 )
	tout << "(";      

      for(int j = 0; j < numCoupledChains; ++j)
	{	  
	  SuccessCounter *ctr = swapInfo[cnt]; 
	  if(i < j )
	    {
	      tout <<  (isFirst ? "" : "," ) << setprecision(1) << 100 * ctr->getRatioOverall() << "%";
	      isFirst = false; 
	    }
	  cnt++; 
	}
      if(i < numCoupledChains - 1 )
	tout << ")";
    }
  
  tout << "\tbaseSwap: " << setprecision(1)<< swapInfo[1]->getRatioInLast100() * 100 << "%,"; 
}


void CoupledChains::switchChainState()
{  
  int numChain = chains.size(); 

  if(numChain == 1)
    return;   

  int chainAId = rand.drawRandInt(numChain),
    chainBId = chainAId; 
  while(chainAId == chainBId)
    chainBId = rand.drawRandInt(numChain); 

  int coupIdA = chains[chainAId].getCouplingId(),
    coupIdB = chains[chainBId].getCouplingId();   

  if(coupIdA > coupIdB)
    swap(coupIdB, coupIdA); 

  double heatA = chains[chainAId].getChainHeat(),
    heatB = chains[chainBId].getChainHeat(); 

  assert(heatA <= 1.f || heatB <= 1.f); 

  double lnlA = chains[chainAId].getTraln().getTr()->likelihood,
    lnlB = chains[chainBId].getTraln().getTr()->likelihood; 

  double lnPrA = chains[chainAId].getPrior().getLnPrior(), 
    lnPrB = chains[chainBId].getPrior().getLnPrior(); 

  double 
    aB = (lnlA + lnPrA ) *  heatB,
    bA = ( lnlB + lnPrB ) *  heatA,
    aA = ( lnlA + lnPrA ) * heatA,
    bB =  ( lnlB + lnPrB ) *  heatB; 

  double accRatio = exp(( aB + bA )  - (aA + bB )); 

  Chain &a = chains[ chainAId],
    &b = chains[ chainBId] ; 

  int r = MIN(coupIdA, coupIdB ); 
  int c = MAX(coupIdA, coupIdB); 
  
  /* do the swap */
  if( rand.drawRandDouble01()  < accRatio)
    {
      a.switchState(b); 
      swapInfo[r * chains.size() + c]->accept(); 
    } 
  else 
    swapInfo[r * chains.size() + c]->reject(); 
}


void CoupledChains::chainInfo()
{
  // print hot chains
  vector<Chain*> sortedChains(chains.size()); 
  for(auto &chain : chains)
    sortedChains[chain.getCouplingId()] = &chain; 

  auto coldChain = sortedChains[0]; 

  Branch fake(0,0,coldChain->getTraln().getTreeLengthExpensive()); 

  tout << "[run: " << runid << "] "  ; 
  tout << "[time " << CLOCK::duration_cast<CLOCK::duration<double> > (CLOCK::system_clock::now()- timeIncrement   ).count()     << "] "; 
  timeIncrement = CLOCK::system_clock::now();   
  tout << "gen: " << coldChain->getGeneration() ; 
  tout <<  "\tTL=" << setprecision(2)<<  fake.getInterpretedLength(coldChain->getTraln()); 
  tout << "\tlnPr(1)=" << coldChain->getPrior().getLnPrior() << "\tlnl(1)=" << setprecision(2)<< coldChain->getTraln().getTr()->likelihood << "\t" ; 

  for(nat i = 1 ; i < chains.size(); ++i)
    {
      auto &chain = sortedChains[i]; 
      double heat = chain->getChainHeat();
      assert(heat < 1.0f); 
      
      tout << "lnl(" << setprecision(2)<< heat << ")=" << setprecision(2)<< chain->getTraln().getTr()->likelihood << "\t" ;  
    }

  printSwapInfo();
  tout << endl; 

  coldChain->printProposalState(tout);

  tout << endl; 
}


void CoupledChains::executePart(int gensToRun)
{  
  /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
  for(auto &c : chains)
    c.resume();

  for(int genCtr = 0; genCtr < gensToRun; genCtr += swapInterval)
    {
      bool timeToPrint = false; 

      for(auto &chain : chains)
	{
	  for(int i = 0; i < swapInterval; ++i)
	    {
	      chain.step(); 	      
	      if(chain.getChainHeat() == 1 
		 && (chain.getGeneration() % samplingFreq)  == samplingFreq - 1
		 ) 
		chain.printSample(topoFile, paramFile); 
		
	      
	      timeToPrint |= 
		chain.getCouplingId() == 0 && printFreq > 0
		&& (chain.getGeneration() % printFreq )  == (printFreq - 1) ;
	    }
	}

      if(timeToPrint)
      	chainInfo(); 

      if(chains.size()  > 1 
	 && tuneHeat
	 && tuneFreq < swapInfo[1]->getRecentlySeen()  )
	{	  
	  tuneTemperature();      	  
	}
      
      switchChainState();
    }

  for(auto &chain : chains)
    chain.suspend();

#ifdef _USE_GOOGLE_PROFILER
  ProfilerFlush();
#endif
}


void CoupledChains::tuneTemperature()
{
  /* naive strategy: tune, s.t. the coldest hot chain swaps
     with the coldest chain in 23.4% of all cases */

  auto c = swapInfo[1]; 
  double deltaT = chains[0].getDeltaT(); 
  deltaT = tuneParameter(  c->getBatch() , c->getRatioInLastInterval(), deltaT, false); 
  c->nextBatch();
  
  // update the chains 
  for(auto& chain : chains)
    chain.setDeltaT(deltaT);   
}

