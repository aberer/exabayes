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


// TODO 
CoupledChains::CoupledChains(randCtr_t seed, int runNum, string workingdir, std::string runname, int numCoupled,  vector<Chain> &_chains   )
  : chains(std::move(_chains))
  , swapInfo(chains.size())
  , heatIncrement(0.1) 
  , rand(seed)
  , runid(runNum) 
  , tuneHeat(false)
  , printFreq(500)
  , swapInterval(1)
  , samplingFreq(100)
  , runname(runname) 
  , workdir(workingdir)
{  
  auto params = chains[0].extractParameters(); 
  nat ctr = 0; 
  for( auto &param : params)
    if(param->getCategory() == Category::BRANCH_LENGTHS)
      ++ctr; 

  for(nat i = 0; i < ctr ; ++i)
    paramId2TopFile.emplace(i, TopologyFile(workingdir, runname, runid,0, i, ctr > 1)); 

  pFile.emplace_back(workdir, runname,runid, 0); 
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
  // , tFile(std::move(rhs.tFile))
  , paramId2TopFile(std::move(rhs.paramId2TopFile))
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
  auto &params = chains[0].extractParameters();

  auto tag =  rand.getKey();

  for(auto &elem : paramId2TopFile)
    elem.second.initialize(traln, tag.v[0]); 
  
  pFile[0].initialize(traln, params, tag.v[0] ); 
  
}


void CoupledChains::seedChains()
{
  for(auto &c : chains)
    c.reseed(rand.generateSeed()); 
  // std::cout << rand  << " after RESEEDING " << std::endl; 
}


void CoupledChains::attemptSwap(ParallelSetup &pl)
{  
  CommFlag flags = CommFlag::PrintStat | CommFlag::Proposals; 

  int numChain = chains.size(); 

  if(numChain == 1)
    return;   

  int cAIndex = rand.drawIntegerOpen(numChain); 
  int cBIndex = rand.drawIntegerOpen(numChain-1) ; 
  if(cBIndex == cAIndex)
    cBIndex = numChain-1; 
  // std::cout << rand  << " attempting swap between indices " << cAIndex << " and " << cBIndex; 
  if( not pl.isMyChain(runid, cAIndex) && not pl.isMyChain(runid,cBIndex))
    {
      rand.drawRandDouble01();	// need to waste one value...  
      // std::cout << " NOT MINE " << std::endl; 
      return;
    }
  // std::cout << std::endl; 

  
  if(not pl.isMyChain(runid, cAIndex))
    std::swap(cAIndex, cBIndex); 

  Chain& a = chains[cAIndex]; 
  Chain& b = chains[cBIndex]; 

  bool mineHasSmallerId = a.getCouplingId() < b.getCouplingId(); 

  auto aSer = a.serializeConditionally( flags ); 
  auto bSer = pl.sendRecvChain(*this, cAIndex, cBIndex, aSer, flags);
  
  if(not pl.isMyChain(runid, cBIndex))
    b.deserializeConditionally(bSer, flags); 

  assert(b.getChainHeat() <= 1. && a.getChainHeat() <= 1.); 

  double aB = (a.getLikelihood() + a.getLnPr()) * b.getChainHeat(),
    bA = (b.getLikelihood() + b.getLnPr()) * a.getChainHeat(),
    aA = (a.getLikelihood() + a.getLnPr())* a.getChainHeat(),
    bB = (b.getLikelihood() + b.getLnPr())* b.getChainHeat();

  double accRatio = exp(( aB + bA )  - (aA + bB )); 

  bool bothAreMine = pl.isMyChain(runid, cAIndex) && pl.isMyChain(runid, cBIndex); 
  nat coupIdA = a.getCouplingId(), 
    coupIdB = b.getCouplingId(); 

  /* do the swap */
  double r = rand.drawRandDouble01(); 
  bool didAccept = r < accRatio;   
  if( didAccept )
    {
      a.deserializeConditionally(bSer,flags); 
      b.deserializeConditionally(aSer, flags); 
    } 
  
  // update swap matrix: if the other chain did not belong to us and
  // our id was greater, do not store the info
  if(bothAreMine || mineHasSmallerId)
    swapInfo.update(coupIdA,coupIdB,didAccept); 
}


void CoupledChains::executePart(nat startGen, nat numGen, ParallelSetup &pl)
{ 

  assert(pl.isMyRun(getRunid())); 

  for(nat i = 0; i < chains.size(); ++i)
    {
      if(pl.isMyChain(runid, i))
	chains[i].resume(true,true);
    }

  // additional sampling, if we are in the very first generation 
  if(startGen == 0)
    {
      for(auto &c : chains)
	if( c.getChainHeat() == 1. && pl.isChainLeader() )
	  c.sample(paramId2TopFile, pFile[0]); 
    }


  for(nat genCtr = startGen; genCtr < startGen + numGen; genCtr += swapInterval)
    {
      for(nat j = 0; j < chains.size(); ++j)
	{
	  auto &chain = chains[j];
	  if(not pl.isMyChain(runid, j))
	    continue; 

	  for(int i = 0; i < swapInterval; ++i)
	    {
	      chain.step(); 	      
	      if(chain.getChainHeat() == 1 && (chain.getGeneration() % samplingFreq)  == 0 ) 
		{
		  for(auto &c : chains)
		    {
		      if( c.getChainHeat() == 1. && pl.isChainLeader() )
			c.sample(paramId2TopFile, pFile[0]); 
		    }
		}
	    }
	}
#ifdef UNSURE  
      assert(0);
#endif
      
      rand.rebase(genCtr); 

      attemptSwap(pl);
    }


  for(nat i = 0; i < chains.size(); ++i)
    {
      auto &chain = chains[i]; 
      if(pl.isMyChain(runid,i))
	chain.suspend();
    }

#ifdef _USE_GOOGLE_PROFILER
  ProfilerFlush();
#endif
}


void CoupledChains::finalizeOutputFiles()   
{
  for(auto &elem : paramId2TopFile)
    elem.second.finalize(); 
  for(auto &p : pFile)
    p.finalize();
}


void CoupledChains::readFromCheckpoint( std::istream &in ) 
{
  rand.readFromCheckpoint(in); 
  swapInfo.readFromCheckpoint(in);
  for(auto &chain : chains)
    chain.readFromCheckpoint(in);
} 


void CoupledChains::writeToCheckpoint( std::ostream &out)   const
{
  rand.writeToCheckpoint(out);
  swapInfo.writeToCheckpoint(out);
  for(auto &chain : chains)
    chain.writeToCheckpoint(out); 
}   


void CoupledChains::regenerateOutputFiles(std::string workdir, std::string prevId) 
{
  nat gen = chains[0].getGeneration();
  for(auto &pF : pFile)
    pF.regenerate(workdir, prevId, gen); 
  for(auto &elem : paramId2TopFile)
    elem.second.regenerate(workdir, prevId, gen); 
} 

