#include <sstream>

#include "CoupledChains.hpp"
#include "Chain.hpp"
#include "adapters.h" 
#include "SuccessCounter.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "treeRead.h"
#include "output.h"
#include "AbstractProposal.hpp"
#include "PriorBelief.hpp"
#include "Category.hpp"

CoupledChains::CoupledChains(int seed, int numCoupled, vector<TreeAln*> trees, int _runid , double _printFreq, double _swapInterval, int _samplingFreq, double heatFactor, string _runname, string workingdir, const PriorBelief &prior, const vector<Category> proposals, int tuneFreq)
  : temperature(heatFactor)
  , rand(seed)
  , runid(_runid) 
  , tuneHeat(false)
  , printFreq(_printFreq)
  , swapInterval(_swapInterval)
  , samplingFreq(_samplingFreq)
  , runname(_runname)
{
  assert((nat)numCoupled == trees.size());

  for(int i = 0; i < numCoupled; ++i)
    {
      vector<Category> pCopy; 
      for(auto c : proposals)
	{
	  Category aCopy(c); 
	  aCopy.copyDeep(c); 
	  pCopy.push_back(aCopy);
	}

      Chain *chain = new Chain(rand.generateSeed(),i, runid, trees[i], prior, pCopy, tuneFreq); 
      chain->setDeltaT(temperature); 
      chains.push_back(chain);
    }

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

  chains[0]->printNexusTreeFileStart(topoFile);
  chains[0]->printParamFileStart(paramFile) ;
}




CoupledChains::~CoupledChains()
{
  // IMPORTANT TODO 
  // for(auto chain : chains) 
  //   exa_free(chain); 

  // for(auto swapI : swapInfo)
  //   exa_free(swapI); 
}


void CoupledChains::printSwapInfo()
{
  if(chains.size( )== 1 )
    return ; 
  
  int numCoupledChains = chains.size(); 
  
  int cnt = 0; 
  for(int i = 0; i < numCoupledChains; ++i)
    {
      if(i < numCoupledChains - 1 )
	tout << "(";


      for(int j = 0; j < numCoupledChains; ++j)
	{
	  SuccessCounter *ctr = swapInfo[cnt]; 
	  if(i < j )
	    tout << setprecision(1) << 100 * ctr->getRatioOverall() << "%";

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

  int coupIdA = chains[chainAId]->getCouplingId(),
    coupIdB = chains[chainBId]->getCouplingId();   

  if(coupIdA > coupIdB)
    swap(coupIdB, coupIdA); 


  /* 
     IMPORTANT TODO

     this currently assumes that we have a non-informative
     prior. Allow for changes!
  */

  double heatA = chains[chainAId]->getChainHeat(),
    heatB = chains[chainBId]->getChainHeat(); 

  assert(heatA <= 1.f || heatB <= 1.f); 

  double lnlA = chains[chainAId]->traln->getTr()->likelihood,
    lnlB = chains[chainBId]->traln->getTr()->likelihood; 

  double 
    aB = lnlA *  heatB,
    bA = lnlB *  heatA,
    aA = lnlA * heatA,
    bB =  lnlB *  heatB; 

  double accRatio = exp(( aB + bA )  - (aA + bB )); 

  Chain *a = chains[ chainAId],
    *b = chains[ chainBId] ; 

  int r = MIN(coupIdA, coupIdB ); 
  int c = MAX(coupIdA, coupIdB); 
  
  /* do the swap */
  if( rand.drawRandDouble01()  < accRatio)
    {
      a->switchState(*b); 
      swapInfo[r * chains.size() + c]->accept(); 
    } 
  else 
    swapInfo[r * chains.size() + c]->reject(); 
}


void CoupledChains::chainInfo()
{
  // find cold chain
  Chain *coldChain = NULL; 
  for(auto chain : chains )
    if(chain->getCouplingId() == 0)
      coldChain = chain; 
  assert(coldChain != NULL); 
  
  tree *tr = coldChain->traln->getTr(); 

  tout << "[run: " << runid << "] [time " << setprecision(2) << gettime()- timeIncrement << "] gen: " << coldChain->getGeneration() <<  "\tTL=" << setprecision(2)<< branchLengthToReal(tr, coldChain->traln->getTreeLength()) << "\tlnPr(1)=" << coldChain->getPrior().getLogProb() << "\tlnl(1)=" << setprecision(2)<< coldChain->traln->getTr()->likelihood << "\t" ; 

  // print hot chains
  vector<Chain*> sortedChains(chains.size()); 
  for(auto chain : chains)
    sortedChains[chain->getCouplingId()] = chain; 
  for(nat i = 1 ; i < chains.size(); ++i)
    {
      Chain *chain = sortedChains[i]; 
      double heat = chain->getChainHeat();
      assert(heat < 1.0f); 
      
      tout << "lnl(" << setprecision(2)<< heat << ")=" << setprecision(2)<< chain->traln->getTr()->likelihood << "\t" ;  
    }

  
  printSwapInfo();


  tout << endl; 

  /* just output how much time has passed since the last increment */
  timeIncrement = gettime(); 	

  for(auto cat : coldChain->getProposalCategories())
    {
      tout << cat.getName() << ":\t";
      for( auto p : cat.getProposals() )
  	tout << p->getName() << ":"  << p->getSCtr() << "\t" ; 
      tout << endl; 
    }

  tout << endl; 
}


void CoupledChains::executePart(int gensToRun)
{  
  /* if we have ample space, then we'll have to use the apply and save functions only at the beginning and end of each run for all chains  */
  for(nat i = 0; i < chains.size(); ++i)
    chains[i]->applyChainStateToTree();

  for(int genCtr = 0; genCtr < gensToRun; genCtr += swapInterval)
    {
      bool timeToPrint = false; 

      for(auto chain : chains)
	{
	  for(int i = 0; i < swapInterval; ++i)
	    {
	      chain->step(); 	      
	      if(chain->getChainHeat() == 1 
		 && (chain->getGeneration() % samplingFreq)  == samplingFreq - 1
		 && isOutputProcess() ) 
		chain->printSample(topoFile, paramFile); 
		
	      
	      timeToPrint |= isOutputProcess() 
		&& chain->getCouplingId() == 0 && printFreq > 0
		&& (chain->getGeneration() % printFreq )  == (printFreq - 1) ;
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

  for(auto chain : chains)
    chain->saveTreeStateToChain();
}


void CoupledChains::tuneTemperature()
{
  /* naive strategy: tune, s.t. the coldest hot chain swaps
     with the coldest chain in 23.4% of all cases */

  auto c = swapInfo[1]; 

  temperature = tuneParameter(  c->getBatch() , c->getRatioInLastInterval(), temperature, FALSE); 
  // tout << "new temperature " << setprecision(3) << temperature<< ". Ratio was "  << c->getRatioInLastInterval() << endl; 
  c->nextBatch();
  
  // update the chains 
  for(auto chain : chains)
    chain->setDeltaT(temperature);   
}

