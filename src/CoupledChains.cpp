#include <sstream>

#include "CoupledChains.hpp"
#include "Chain.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "treeRead.h"
#include "AbstractProposal.hpp"
#include "PriorBelief.hpp"

CoupledChains::CoupledChains(int seed, int runNum, const BlockRunParameters &params, vector<TreeAlnPtr > trees, string workingdir, const vector<ProposalPtr> &proposals, const vector<RandomVariablePtr> &vars, LikelihoodEvaluatorPtr eval)
  : temperature(params.getHeatFactor())
  , rand(seed)
  , runid(runNum) 
  , tuneHeat(params.getTuneHeat())
  , printFreq(params.getPrintFreq())
  , swapInterval(params.getSwapInterval())
  , samplingFreq(params.getSamplingFreq())
  , runname(params.getRunId())
{
  int numCoupled = params.getNumCoupledChains();  
  assert((nat)numCoupled == trees.size());

  for(int i = 0; i < numCoupled; ++i)
    {

      Chain *chain = new Chain(rand.generateSeed(),i, runid, trees[i],  proposals, params.getTuneFreq(), vars, eval); 
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

  int coupIdA = chains[chainAId]->getCouplingId(),
    coupIdB = chains[chainBId]->getCouplingId();   

  if(coupIdA > coupIdB)
    swap(coupIdB, coupIdA); 

  double heatA = chains[chainAId]->getChainHeat(),
    heatB = chains[chainBId]->getChainHeat(); 

  assert(heatA <= 1.f || heatB <= 1.f); 

  double lnlA = chains[chainAId]->getTraln().getTr()->likelihood,
    lnlB = chains[chainBId]->getTraln().getTr()->likelihood; 

  double lnPrA = chains[chainAId]->getPrior().getLnPrior(), 
    lnPrB = chains[chainBId]->getPrior().getLnPrior(); 

  double 
    aB = (lnlA + lnPrA ) *  heatB,
    bA = ( lnlB + lnPrB ) *  heatA,
    aA = ( lnlA + lnPrA ) * heatA,
    bB =  ( lnlB + lnPrB ) *  heatB; 

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
  // print hot chains
  vector<Chain*> sortedChains(chains.size()); 
  for(auto chain : chains)
    sortedChains[chain->getCouplingId()] = chain; 

  const Chain& coldChain = *(sortedChains[0]); 

  Branch fake(0,0,coldChain.getTraln().getTreeLengthExpensive()); 

  tout << "[run: " << runid << "] [time " << setprecision(2) << gettime()- timeIncrement << "] gen: " << coldChain.getGeneration() 
       <<  "\tTL=" << setprecision(2)<<  fake.getInterpretedLength(coldChain.getTraln())
       << "\tlnPr(1)=" << coldChain.getPrior().getLnPrior() << "\tlnl(1)=" << setprecision(2)<< coldChain.getTraln().getTr()->likelihood << "\t" ; 

  for(nat i = 1 ; i < chains.size(); ++i)
    {
      Chain *chain = sortedChains[i]; 
      double heat = chain->getChainHeat();
      assert(heat < 1.0f); 
      
      tout << "lnl(" << setprecision(2)<< heat << ")=" << setprecision(2)<< chain->getTraln().getTr()->likelihood << "\t" ;  
    }

  printSwapInfo();
  tout << endl; 
  timeIncrement = gettime(); 	

  map<Category, vector< AbstractProposal* > > sortedProposals; 
  auto& proposals = coldChain.getProposals();
  for(auto& p : proposals)
    sortedProposals[p->getCategory()].push_back(p.get()) ; 

  for(auto &n : CategoryFuns::getAllCategories())
    {       
      Category cat = n; 

      bool isThere = false; 
      for(auto &p : proposals)
	isThere |= p->getCategory() == cat ; 
      
      if(isThere)
	{
	  tout << CategoryFuns::getLongName(n) << ":\t";
	  for(auto &p : proposals)
	    {
	      if(p->getCategory() == cat)
		{
		  p->printNamePartitions(tout); 
		  tout  << ":"  << p->getSCtr() << "\t" ; 	      
		}
	    }
	  tout << endl; 
	}
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
		 ) 
		chain->printSample(topoFile, paramFile); 
		
	      
	      timeToPrint |= 
		chain->getCouplingId() == 0 && printFreq > 0
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

#ifdef _USE_GOOGLE_PROFILER
  ProfilerFlush();
#endif
}


void CoupledChains::tuneTemperature()
{
  /* naive strategy: tune, s.t. the coldest hot chain swaps
     with the coldest chain in 23.4% of all cases */

  auto c = swapInfo[1]; 

  temperature = tuneParameter(  c->getBatch() , c->getRatioInLastInterval(), temperature, false); 
  // tout << "new temperature " << setprecision(3) << temperature<< ". Ratio was "  << c->getRatioInLastInterval() << endl; 
  c->nextBatch();
  
  // update the chains 
  for(auto chain : chains)
    chain->setDeltaT(temperature);   
}

