#include <sstream>

#include "PendingSwap.hpp"
#include "CoupledChains.hpp"   
#include "Chain.hpp"
#include "GlobalVariables.hpp"
#include "tune.h"
#include "proposals/AbstractProposal.hpp"
#include "PriorBelief.hpp"
#include "time.hpp"
#include "ParallelSetup.hpp"

#define USE_NONBLOCKING_COMM

// #define VERBOSE 

CoupledChains::CoupledChains(Randomness randI, int runNum, string workingdir, std::string runname, int numCoupled,  std::vector<Chain> &chains   )
  : _chains(std::move(chains))
  , _swapInfo(_chains.size())
  , _heatIncrement(0.1) 
  , _rand(randI)
  , _runid(runNum) 
  , _samplingFreq(100)
  , _runname(runname) 
  , _workdir(workingdir)
  , _numSwapsPerGen{1.}
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


PendingSwap CoupledChains::prepareSwap(ParallelSetup &pl, const SwapElem& theSwap)
{
  auto flags = CommFlag::PRINT_STAT | CommFlag::PROPOSALS; 

  auto myId = theSwap.getMyId(pl,_runid); 

  auto &&oss = std::ostringstream{};
  _chains[myId].serializeConditionally(oss,flags);

  auto mySer = std::vector<char>{}; 
  auto str = oss.str();
  mySer.insert(end(mySer), begin(str), end(str));

  auto result = PendingSwap{theSwap}; 
  result.initialize(pl, mySer, _runid, _chains.size());
  return result; 
}


bool CoupledChains::doSwap(ParallelSetup &pl, const SwapElem &theSwap )
{  
  int numChain = _chains.size(); 
  assert(numChain > 1 ); 

  auto flags = CommFlag::PRINT_STAT | CommFlag::PROPOSALS; 

  int cAIndex = theSwap.getOne() ; 
  int cBIndex = theSwap.getOther()  ; 

  assert(  pl.isMyChain(_runid, cAIndex) ||    pl.isMyChain(_runid,cBIndex)) ; 

  if( not pl.isMyChain(_runid, cAIndex) )
    std::swap(cAIndex, cBIndex); 

  auto& a = _chains[cAIndex];
  auto& b = _chains[cBIndex];

  bool mineHasSmallerId = a.getCouplingId() < b.getCouplingId();
  bool bothAreMine = pl.isMyChain(_runid, cAIndex) && pl.isMyChain(_runid, cBIndex);
  
  auto aSer = std::string{}; 
  auto bSer = std::string{}; 

  if(not bothAreMine)
    {
      auto &&ass = std::ostringstream{}; 
      a.serializeConditionally( ass, flags ); 
      bSer = pl.sendRecvChain(*this, cAIndex, cBIndex, ass.str(), flags);
      auto &&iss = std::istringstream{bSer};
      b.deserializeConditionally(iss, flags); 
    }

  assert(a.getGeneration() == b.getGeneration()); 
  assert(b.getChainHeat() <= 1. && a.getChainHeat() <= 1.); 

  double aB = (a.getLikelihood() + a.getLnPr()) * b.getChainHeat(),
    bA = (b.getLikelihood() + b.getLnPr()) * a.getChainHeat(),
    aA = (a.getLikelihood() + a.getLnPr()) * a.getChainHeat(),
    bB = (b.getLikelihood() + b.getLnPr()) * b.getChainHeat();

  
  double accRatio = min(exp(( aB + bA )  - (aA + bB )),1.0); 

  nat coupIdA = a.getCouplingId(), 
    coupIdB = b.getCouplingId(); 

  /* do the swap */
  double r = theSwap.getR();
  bool didAccept = r < accRatio;   
  if( didAccept )
    {
      swapHeatAndProposals(a,b);
    }
  
  // update swap matrix: if the other chain did not belong to us and
  // our id was greater, do not store the info
  if(bothAreMine || mineHasSmallerId)
    _swapInfo.update(coupIdA,coupIdB,didAccept); 

  return didAccept; 
}


bool CoupledChains::doLocalSwap(ParallelSetup &pl, const SwapElem &theSwap )
{  
  int numChain = _chains.size(); 
  assert(numChain > 1 ); 

  int cAIndex = theSwap.getOne() ; 
  int cBIndex = theSwap.getOther()  ; 

  assert(  pl.isMyChain(_runid, cAIndex) ||    pl.isMyChain(_runid,cBIndex)) ; 

  if( not pl.isMyChain(_runid, cAIndex) )
    std::swap(cAIndex, cBIndex); 

  auto& a = _chains[cAIndex];
  auto& b = _chains[cBIndex];

  bool mineHasSmallerId = a.getCouplingId() < b.getCouplingId();
  bool bothAreMine = pl.isMyChain(_runid, cAIndex) && pl.isMyChain(_runid, cBIndex);

  assert(a.getGeneration() == b.getGeneration()); 
  assert(b.getChainHeat() <= 1. && a.getChainHeat() <= 1.); 

  double aB = (a.getLikelihood() + a.getLnPr()) * b.getChainHeat(),
    bA = (b.getLikelihood() + b.getLnPr()) * a.getChainHeat(),
    aA = (a.getLikelihood() + a.getLnPr()) * a.getChainHeat(),
    bB = (b.getLikelihood() + b.getLnPr()) * b.getChainHeat();

  
  double accRatio = min(exp(( aB + bA )  - (aA + bB )),1.0); 

  nat coupIdA = a.getCouplingId(), 
    coupIdB = b.getCouplingId(); 

  /* do the swap */
  double r = theSwap.getR();
  bool didAccept = r < accRatio;   
  if( didAccept )
    {
      swapHeatAndProposals(a,b);
    }
  
  // update swap matrix: if the other chain did not belong to us and
  // our id was greater, do not store the info
  if(bothAreMine || mineHasSmallerId)
    _swapInfo.update(coupIdA,coupIdB,didAccept); 

  return didAccept; 
}


static bool allFinished(const std::vector<Chain> &chains, nat gen, const ParallelSetup &pl,  nat runid)
{
  bool allFinished = true; 
  for(nat i = 0; i < chains.size() ; ++i)
    allFinished &= not pl.isMyChain(runid, i) || (chains[i].getGeneration() == gen); 
  return allFinished; 
}


std::list<SwapElem>
CoupledChains::generateSwapsForBatch(nat startGen, nat numGen) 
{
  auto allSwaps = std::list<SwapElem>{}; 
  if( _chains.size() > 1  )
    {
      for(nat i = startGen; i < startGen + numGen; ++i)
	{
	  nat gen = i + 1; 
	  _rand.rebaseForGeneration(gen); 
	  
	  double n = 2 * _chains.size() * _numSwapsPerGen; 
	  if(n == 0)
	    n = _chains.size(); 
	  double p = _numSwapsPerGen / n; 

	  nat numSwaps = _rand.drawBinomial(p,n); 
	  
	  for(nat j = 0; j < numSwaps; ++j)
	    {
	      auto cAIndex = _rand.drawIntegerOpen(_chains.size()) ; 
	      auto cBIndex = _rand.drawIntegerOpen(_chains.size()-1); 
	      if(cAIndex == cBIndex)
		cBIndex = _chains.size()-1; 

	      if(cAIndex > cBIndex)
		std::swap(cAIndex, cBIndex); 

	      double r = _rand.drawRandDouble01();
	      allSwaps.emplace_back(gen,cAIndex, cBIndex,r); 
	    }
	}
    }
  return allSwaps; 
}


void CoupledChains::doStep(nat id, ParallelSetup &pl)
{
  auto &chain = _chains[id]; 
  chain.step();
#ifdef VERBOSE
  std::cout << "STEP " << id << "\t" << chain.getGeneration() << std::endl; 
#endif
  if(chain.getChainHeat() == 1.
     && chain.getGeneration() % _samplingFreq == 0 
     && pl.isChainLeader())
    chain.sample(_paramId2TopFile, _pFile[0]);
}



bool CoupledChains::allMyChainsAreBlocked( const std::vector<bool> &isBlocked, const ParallelSetup& pl ) const 
{
  auto result = true;  
  for(nat i = 0; i < isBlocked.size();  ++i)
    {
      if(pl.isMyChain(_runid,i))
	result &= isBlocked.at(i); 
    }
  return result; 
}


static void block(std::vector<bool> &isBlocked, nat num)
{
  isBlocked[num] = true; 
}

static void unblock(std::vector<bool> &isBlocked, nat num)
{
  isBlocked[num] = false; 
}




void CoupledChains::executePartNew(nat startGen, nat numGen, ParallelSetup& pl)
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
      nat ctr = 0;  
      for(auto &c : _chains)
  	if( c.getChainHeat() == 1. && pl.isChainLeader() && pl.isMyChain(_runid, ctr) )
	  {
	    c.sample(_paramId2TopFile, _pFile[0]); 
	    ++ctr; 
	  }
    }

  auto allSwaps = generateSwapsForBatch(startGen, numGen);

  for(nat i = 0; i < _chains.size() ; ++i)
    {
      if(pl.isMyChain(_runid, i))
	assert(_chains[i].getGeneration() == startGen ); 
    }

  auto pendingSwaps = std::list<PendingSwap>{}; 
  auto swapsToBeDeleted = std::list<PendingSwap>{}; 
#ifdef USE_NONBLOCKING_COMM
  auto flags = CommFlag::PRINT_STAT | CommFlag::PROPOSALS; 
  while(not allSwaps.empty() ||  not pendingSwaps.empty())
    {
      auto isBlocked = std::vector<bool>( _chains.size(),false ); 
      
      for(auto &elem : pendingSwaps)
	{
	  block(isBlocked, elem.getSwap().getOne()) ;
	  block(isBlocked, elem.getSwap().getOther()); 
	}


#ifdef VERBOSE
      std::cout << "PENDING:\t" ; 
      for(auto &elem : pendingSwaps)
      	std::cout << "{" << elem.getSwap() << "}\t" ; 
      std::cout << std::endl; 

      std::cout << "BLOCKED:\t"; 
      for(auto elem : isBlocked)
      	std::cout << ( elem ? "true" : "false")  << ","; 
      std::cout << std::endl; 
#endif

      // check, if a request completed 
      {
	auto psIter = begin(pendingSwaps); 
	while(psIter != end(pendingSwaps))
	  {
	    if(psIter->allHaveReceived(pl))
	      {
		auto serializedRemote = psIter->getRemoteData();
		
		auto swap = psIter->getSwap(); 
		auto remoteId = swap.getRemoteId(pl,_runid);
		auto myId = swap.getMyId(pl,_runid);

#ifdef VERBOSE
		std::cout << "swap finished: " << swap << std::endl;
#endif

		auto &&iss = std::istringstream{serializedRemote}; 
		_chains[remoteId].deserializeConditionally(iss, flags);

		assert(_chains[remoteId].getGeneration() == swap.getGen() ); 
		assert(_chains[myId].getGeneration() == swap.getGen()); 

		doLocalSwap(pl, swap);

		unblock(isBlocked, myId);

		auto &&elem = std::move(*psIter); 
		swapsToBeDeleted.insert(end(swapsToBeDeleted), std::move(elem));

		psIter = pendingSwaps.erase(psIter);
	      }
	    else 
	      ++psIter; 
	  }
      }

      // check, what can be deleted
      {
	auto psIter = begin(swapsToBeDeleted); 
	while(psIter != end(swapsToBeDeleted))
	  {
	    if(psIter->isFinished())
	      psIter = swapsToBeDeleted.erase(psIter); 
	    else 
	      ++psIter; 
	  }
      }

      auto swapIter = begin(allSwaps); 
      while(swapIter != end(allSwaps) && not allMyChainsAreBlocked(isBlocked,pl))
	{
	  auto swap = *swapIter; 

	  auto one = swap.getOne(); 
	  auto other = swap.getOther();
	  auto gen = swap.getGen();

	  auto bothAreMine = swap.bothAreMine(pl,_runid); 
	  auto oneIsMine = swap.onlyOneIsMine(pl, _runid);
	  
	  if(not bothAreMine && not oneIsMine )
	    {
	      swapIter = allSwaps.erase(swapIter); 
	      continue; 
	    }
	  else if(oneIsMine)
	    {
	      nat myId = swap.getMyId(pl, _runid);
	      while(  not isBlocked[myId] &&  _chains[myId].getGeneration() < gen)
		doStep(myId, pl); 

	      if(not isBlocked[myId] && gen == _chains[myId].getGeneration())
		{
		  pendingSwaps.insert(end(pendingSwaps), std::move(prepareSwap(pl, swap))); 
#ifdef VERBOSE
		  std::cout << "QUEUED " << swap << std::endl; 
#endif
		  swapIter = allSwaps.erase(swapIter);
		}
	      else 
		{
		  ++swapIter; 
		}
	      
	      block(isBlocked, myId); 
	    }
	  else 
	    {
	      while( not isBlocked[one] && _chains[one].getGeneration() < gen)
		doStep(one,pl); 
	      while( not isBlocked[other] && _chains[other].getGeneration() < gen) 
		doStep(other,pl); 

	      if(not isBlocked[one] && not isBlocked[other]  
		 && gen == _chains[one].getGeneration() && gen == _chains[other].getGeneration())
		{
		  doSwap(pl, swap);
		  swapIter = allSwaps.erase(swapIter);
#ifdef VERBOSE
		  std::cout << "LOCAL " << swap << std::endl; 
#endif
		}
	      else 
		{
		  block(isBlocked, one); 
		  block(isBlocked, other); 

		  ++swapIter; 
		}
	    }
	}

#ifdef VERBOSE      
      std::cout << std::endl; 
      sleep(1); 
#endif
    }
#else 
  // old version 
  while(not allSwaps.empty() || not pendingSwaps.empty())
    {
      auto swap = allSwaps.front(); 
      allSwaps.erase(begin(allSwaps));

      auto one = swap.getOne(); 
      auto other = swap.getOther();
      auto gen = swap.getGen();
      
      auto bothAreMine = pl.isMyChain(_runid, one) && pl.isMyChain(_runid,other); 
      auto oneIsMine = ( not pl.isMyChain(_runid, one) ) != ( not  pl.isMyChain(_runid,other)) ; 

      if(not bothAreMine && not oneIsMine )
	continue; 

      if(oneIsMine)
	{
	  nat myId = swap.getMyId(pl,_runid); 
	  while( _chains[myId].getGeneration() < gen)
	    doStep(myId, pl); 
	}
      else 
	{
	  while( _chains[one].getGeneration() < gen)
	    doStep(one,pl); 
	  while( _chains[other].getGeneration() < gen) 
	    doStep(other,pl); 
	}

      doSwap(pl, swap);
    }
#endif

  // execute remaining steps
  for(nat i = 0; i < _chains.size() ; ++i)
    {
      if(pl.isMyChain(_runid, i))
	{
	  while(_chains[i].getGeneration() < startGen + numGen)
	    doStep(i,pl);
	}
    }

  for(nat i = 0; i < _chains.size(); ++i)
    {
      auto &chain = _chains[i]; 
      
      if(pl.isMyChain(_runid,i))
  	{
  	  assert(chain.getGeneration() == startGen + numGen); 
  	  chain.suspend();
  	}
    }
}


void CoupledChains::finalizeOutputFiles(const ParallelSetup &pl) 
{
  if(pl.isRunLeader())
    {
      for(auto &elem : _paramId2TopFile)
	elem.second.finalize(); 
      for(auto &p : _pFile)
	p.finalize();
    }
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
