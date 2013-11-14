#include <sstream>

#include "CoupledChains.hpp"   
#include "Chain.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "proposals/AbstractProposal.hpp"
#include "PriorBelief.hpp"
#include "time.hpp"
#include "ParallelSetup.hpp"


CoupledChains::CoupledChains(Randomness randI, int runNum, string workingdir, std::string runname, int numCoupled,  std::vector<Chain> &chains   )
  : _chains(std::move(chains))
  , _swapInfo(_chains.size())
  , _heatIncrement(0.1) 
  , _rand(randI)
  , _runid(runNum) 
  , _samplingFreq(100)
  , _runname(runname) 
  , _workdir(workingdir)
  , _numSwapsPerGen(1.)
{  
  auto params = _chains[0].extractParameters(); 

  
  AbstractParameter* topoParamUnfixed = nullptr; 
  auto blParamsUnfixed = std::vector<AbstractParameter*>{}; 
  for( auto &param : params)
    {
      if(param->getCategory() == Category::BRANCH_LENGTHS && param->getPrior()->needsIntegration()  )
	blParamsUnfixed.push_back(param); 
      else if(param->getCategory() == Category::TOPOLOGY && param->getPrior()->needsIntegration() ) 
	topoParamUnfixed = param;
    }

  
  if(blParamsUnfixed.size() > 0)
    {
      nat ctr = 0; 
      for(auto &param : blParamsUnfixed)
	{
	  _paramId2TopFile.insert(std::make_pair(param->getId(), TopologyFile(workingdir, _runname, _runid,0, ctr, blParamsUnfixed.size() > 1))); 
	  ++ctr; 
	}
    }
  else if(topoParamUnfixed != nullptr)
    _paramId2TopFile.insert(std::make_pair(topoParamUnfixed->getId(), TopologyFile(workingdir, _runname, _runid,0, 0, false))); 

  _pFile.emplace_back(_workdir, _runname,_runid, 0); 
}


CoupledChains::CoupledChains(CoupledChains&& rhs)   
  : _chains(std::move(rhs._chains))
  , _swapInfo(std::move(rhs._swapInfo))    
  , _heatIncrement(rhs._heatIncrement)
  , _rand(std::move(rhs._rand))
  , _runid(rhs._runid)
  , _samplingFreq(rhs._samplingFreq)
  , _runname(std::move(rhs._runname))
  , _workdir(std::move(_workdir))
  , _paramId2TopFile(std::move(rhs._paramId2TopFile))
  , _pFile(std::move(rhs._pFile))
  , _numSwapsPerGen(rhs._numSwapsPerGen)
{
  
}

CoupledChains& CoupledChains::operator=(CoupledChains rhs)
{
  std::swap(*this, rhs); 
  return *this; 
} 


void CoupledChains::initializeOutputFiles(bool isDryRun)  
{  
  // TODO sampling file for every chain possibly 
  auto &traln = _chains[0].getTralnHandle(); 
  auto params = _chains[0].extractParameters();

  auto tag =  _rand.getKey();

  for(auto &elem : _paramId2TopFile)
    elem.second.initialize(traln, tag.v[0], isDryRun); 
  
  _pFile[0].initialize(traln, params, tag.v[0] ,isDryRun); 
  
}


void CoupledChains::attemptSwap(ParallelSetup &pl)
{  
  auto flags = CommFlag::PrintStat | CommFlag::Proposals; 

  int numChain = _chains.size(); 

  if(numChain == 1)
    return;   

  int cAIndex = _rand.drawIntegerOpen(numChain); 
  int cBIndex = _rand.drawIntegerOpen(numChain-1) ; 
  if(cBIndex == cAIndex)
    cBIndex = numChain-1; 
  if( not pl.isMyChain(_runid, cAIndex) && not pl.isMyChain(_runid,cBIndex))
    {
      _rand.drawRandDouble01();	// need to waste one value...  
      return;
    }

  if(not pl.isMyChain(_runid, cAIndex))
    std::swap(cAIndex, cBIndex); 

  auto& a = _chains[cAIndex]; 
  auto& b = _chains[cBIndex]; 

  bool mineHasSmallerId = a.getCouplingId() < b.getCouplingId(); 
  bool bothAreMine = pl.isMyChain(_runid, cAIndex) && pl.isMyChain(_runid, cBIndex); 
  
  auto aSer = std::string{}; 
  auto bSer = std::string{}; 

  if(not bothAreMine)
    {
      aSer = a.serializeConditionally( flags ); 
      bSer = pl.sendRecvChain(*this, cAIndex, cBIndex, aSer, flags);
  
      if(not pl.isMyChain(_runid, cBIndex))
	b.deserializeConditionally(bSer, flags); 
    }

  assert(b.getChainHeat() <= 1. && a.getChainHeat() <= 1.); 

  double aB = (a.getLikelihood() + a.getLnPr()) * b.getChainHeat(),
    bA = (b.getLikelihood() + b.getLnPr()) * a.getChainHeat(),
    aA = (a.getLikelihood() + a.getLnPr()) * a.getChainHeat(),
    bB = (b.getLikelihood() + b.getLnPr()) * b.getChainHeat();

  double accRatio = min(exp(( aB + bA )  - (aA + bB )),1.0); 

  nat coupIdA = a.getCouplingId(), 
    coupIdB = b.getCouplingId(); 

  /* do the swap */
  double r = _rand.drawRandDouble01(); 
  bool didAccept = r < accRatio;   
  if( didAccept )
    {
      if( bothAreMine)
      	{
	  // tout << a.getCouplingId()  << "," << b.getCouplingId() << " => " ; 
	  // std::swap(a,b);
	  // tout << "\t" << a.getCouplingId()  << "," << b.getCouplingId() << std::endl; 
      	  // std::swap(_chains[cAIndex], _chains[cBIndex]); 
	  
	  swapHeatAndProposals(a,b);

      	}
      else 
      	{
	  double lnlA = a.getLikelihood(), 
	    lnPrA = a.getLnPr() ,
	    lnlB = b.getLikelihood(), 
	    lnPrB = b.getLnPr(); 

	  a.deserializeConditionally(bSer,flags); 
	  b.deserializeConditionally(aSer, flags); 

	  // BAAAAD 
	  a.setLikelihood(lnlA); 
	  b.setLikelihood(lnlB); 
	  a.setLnPr(lnPrA); 
	  b.setLnPr(lnPrB); 
	}
    } 

  // update swap matrix: if the other chain did not belong to us and
  // our id was greater, do not store the info
  if(bothAreMine || mineHasSmallerId)
    _swapInfo.update(coupIdA,coupIdB,didAccept); 
}


void CoupledChains::executePart(nat startGen, nat numGen, ParallelSetup &pl)
{ 

  assert(pl.isMyRun(getRunid())); 

  for(nat i = 0; i < _chains.size(); ++i)
    {
      if(pl.isMyChain(_runid, i))
	_chains[i].resume();
    }

  // additional sampling, if we are in the very first generation 
  if(startGen == 0)
    {
      for(auto &c : _chains)
	if( c.getChainHeat() == 1. && pl.isChainLeader() )
	  c.sample(_paramId2TopFile, _pFile[0]); 
    }

  nat endGen = startGen + numGen; 
  nat genCtr = startGen; 
  while(genCtr < endGen)
    {
      nat gensInARow = 0;
      nat numSwaps = 0;  

      while(genCtr + gensInARow < endGen && numSwaps == 0)
	{
	  ++gensInARow;
	  _rand.rebaseForGeneration(genCtr + gensInARow);
	  numSwaps = _rand.drawBinomial(0.5, 2 * _numSwapsPerGen); // HARD-CODED, make parameter
	}
      _rand.rebaseForGeneration(genCtr + gensInARow); 

      for(nat i = 0; i < _chains.size() ;++i)
	{
	  if(not pl.isMyChain(_runid, i))
	    continue; 

	  auto &chain = _chains[i]; 

	  for(nat j = 0; j < gensInARow; ++j)
	    {
	      chain.step();
	      if(chain.getChainHeat() == 1. && (chain.getGeneration() % _samplingFreq)  == 0 ) 
		{
		  for(auto &c : _chains)
		    {
		      if( c.getChainHeat() == 1. && pl.isChainLeader() )
			c.sample(_paramId2TopFile, _pFile[0]); 
		    }
		}
	    }
	}

      genCtr += gensInARow; 
      assert(genCtr == _rand.getGeneration() ); 

      // tout << "swaps=" << numSwaps << " at " << genCtr << std::endl ; 
      
      for(nat i = 0; i < numSwaps; ++i)
	attemptSwap(pl);

      // check 
      for(nat i = 0; i < _chains.size() ; ++i)
	{
	  if(pl.isMyChain(_runid,i))
	    assert(_chains[i].getGeneration() == int(genCtr)); 
	}
    }

  for(nat i = 0; i < _chains.size(); ++i)
    {
      auto &chain = _chains[i]; 
      if(pl.isMyChain(_runid,i))
	chain.suspend();
    }

#ifdef _USE_GOOGLE_PROFILER
  ProfilerFlush();
#endif
}


void CoupledChains::finalizeOutputFiles()   
{
  for(auto &elem : _paramId2TopFile)
    elem.second.finalize(); 
  for(auto &p : _pFile)
    p.finalize();
}


void CoupledChains::deserialize( std::istream &in ) 
{
  _rand.deserialize(in); 
  _swapInfo.deserialize(in);
  for(auto &chain : _chains)
    chain.deserialize(in);
} 


void CoupledChains::serialize( std::ostream &out)   const
{
  _rand.serialize(out);
  _swapInfo.serialize(out);
  for(auto &chain : _chains)
    chain.serialize(out); 
}   


void CoupledChains::regenerateOutputFiles(std::string workdir, std::string prevId) 
{
  nat gen = _chains[0].getGeneration();
  for(auto &pF : _pFile)
    pF.regenerate(workdir, prevId, gen); 
  for(auto &elem : _paramId2TopFile)
    elem.second.regenerate(workdir, prevId, gen); 
} 



std::vector<std::string> CoupledChains::getAllFileNames() const 
{
  auto result = std::vector<std::string>{}; 
  for(auto &elem : _paramId2TopFile)
    result.push_back(elem.second.getFileName());
  for(auto &elem : _pFile)
    result.push_back(elem.getFileName()); 
  
  return result; 
} 
