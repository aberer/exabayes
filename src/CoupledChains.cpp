#include <sstream>

#include "CoupledChains.hpp"   
#include "Chain.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "treeRead.h"
#include "AbstractProposal.hpp"
#include "PriorBelief.hpp"
#include "time.hpp"
#include "ParallelSetup.hpp"

CoupledChains::CoupledChains(randCtr_t seed, int runNum, string workingdir, int numCoupled,  vector<Chain> &_chains   )
  : chains(std::move(_chains))
  , heatIncrement(0.1) 
  , rand(seed)
  , runid(runNum) 
  , tuneHeat(false)
  , printFreq(500)
  , swapInterval(1)
  , samplingFreq(100)
  , runname("standardId") 
  , workdir(workingdir)
{
  // swap info matrix 
  for(int i = 0; i < numCoupled * numCoupled ; ++i)
    swapInfo.push_back(new SuccessCounter()); 
}


CoupledChains::CoupledChains(CoupledChains&& rhs)   
  : chains(std::move(rhs.chains))
  , swapInfo(std::move(rhs.swapInfo))    
  , heatIncrement(rhs.heatIncrement)
  , rand(std::move(rhs.rand))
  , runid(rhs.runid)
  , tuneHeat(rhs.tuneHeat)
  , printFreq(rhs.printFreq)
  , swapInterval(rhs.swapInterval)
  , samplingFreq(rhs.samplingFreq)
  , runname(std::move(rhs.runname))
  , workdir(std::move(workdir))
  , tFile(std::move(rhs.tFile))
  , pFile(std::move(rhs.pFile))
{
  
}

CoupledChains& CoupledChains::operator=(CoupledChains rhs)
{
  std::swap(*this, rhs); 
  return *this; 
} 


void CoupledChains::initializeOutputFiles()  
{

  // TODO sampling file for every chain possibly 
  auto &traln = chains[0].getTraln(); 
  auto &params = chains[0].extractVariables();

  tFile.emplace_back(workdir, runname, runid, 0); 
  pFile.emplace_back(workdir, runname,runid, 0); 
  
  nat id = rand(); 		// meh 
  tFile[0].initialize(traln, id ); 
  pFile[0].initialize(traln, params, id ); 
}


void CoupledChains::seedChains()
{
  for(auto &c : chains)
    c.reseed(rand.generateSeed()); 
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


void CoupledChains::executePart(int gensToRun, const ParallelSetup &pl)
{  
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
	      if(chain.getChainHeat() == 1 && (chain.getGeneration() % samplingFreq)  == samplingFreq - 1 ) 
		{
		  for(auto &c : chains)
		    {
		      if( c.getChainHeat() == 1. && pl.isReportingProcess() )
			c.sample(tFile[0], pFile[0]); 
		    }
		}
		
	      timeToPrint |= 
		chain.getCouplingId() == 0 && printFreq > 0
		&& (chain.getGeneration() % printFreq )  == (printFreq - 1) ;
	    }
	}

#ifdef PRINT_MUCH
      if(timeToPrint)
      	chainInfo(); 
#endif

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


void CoupledChains::finalizeOutputFiles() const 
{
  for(auto &t : tFile)
    t.finalize();
  for(auto &p : pFile)
    p.finalize();
}
