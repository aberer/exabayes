#include "PendingSwap.hpp"

#if HAVE_PLL == 0 

#include "ParallelSetup.hpp"

PendingSwap::PendingSwap(SwapElem swap)
  : _swap{swap}
{
}


PendingSwap::PendingSwap(PendingSwap&& rhs)
  : _sentReqs{std::move(rhs._sentReqs)}
  , _recvReq{std::move(rhs._recvReq)}
  , _swap{std::move(rhs._swap)}
{
}


void swap(PendingSwap& lhs, PendingSwap& rhs)
{
  using std::swap; 
  swap(lhs._sentReqs, rhs._sentReqs);
  swap(lhs._recvReq, rhs._recvReq);
  swap(lhs._swap, rhs._swap);
}


std::string PendingSwap::getRemoteData()  
{
  assert(_recvReq.isServed()); 
  auto remoteData = _recvReq.getArray(); 
  auto result = std::string{begin(remoteData), end(remoteData)}; 
  return result; 
} 


int PendingSwap::createTag( nat numChains ) const 
{
  auto chainA = _swap.getOne(); 
  auto chainB = _swap.getOther(); 
  return std::min(chainA, chainB) * numChains + std::max(chainA, chainB); 
}


void PendingSwap::initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid, nat numChains) 
{
  int tag = createTag( numChains ); 

  nat myId = _swap.getMyId(pl,runid); 
  nat remoteId = _swap.getRemoteId(pl,runid);

  nat ourRunBatch = runid % pl.getRunsParallel();
  nat myBatch = myId % pl.getChainsParallel();
  nat remoteBatch = remoteId % pl.getChainsParallel(); 
  
  nat numInMyBatch = pl.getMesh().getNumRanksInDim(ourRunBatch, myBatch); 
  nat numInRemoteBatch = pl.getMesh().getNumRanksInDim(ourRunBatch, remoteBatch); 

  auto myCoords = pl.getMyCoordinates(); 

  // determine from who i receive 
  nat srcRankInRun = pl.getRankInRun( { {ourRunBatch, remoteBatch, myCoords[2] % numInRemoteBatch }  } ) ; 

  nat myRankInRun =  pl.getRankInRun( myCoords ); 

  // auto &&ss = std::stringstream{}; 
  // ss << "my rank in run is " << myRankInRun << " and in the runComm is " << pl.getRunComm().getRank() << "\tcoords=" << myCoords[0] << "," << myCoords[1] << "," << myCoords[2]; 
  // pl.blockingPrint(pl.getRunComm(), ss.str());
  // pl.globalBarrier();

  assert( int(myRankInRun) == pl.getRunComm().getRank() ); 

  _recvReq.initialize(false, srcRankInRun, tag, std::vector<char>(myChainSer.size()), pl.getRunComm());

  // send stuff 
  for(nat i = 0 ; i < numInRemoteBatch; ++i)
    {
      if( i % numInMyBatch == myCoords[2])
	{
	  nat hisRank = pl.getRankInRun( {{ ourRunBatch, remoteBatch, i } }); 

	  auto &&req =  CommRequest(); 
	  _sentReqs.push_back(std::move(req));
	  _sentReqs.back().initialize(true, hisRank, tag, myChainSer, pl.getRunComm());
	}
    }
}


bool PendingSwap::allHaveReceived(ParallelSetup& pl)  
{
  bool hasReceived = _recvReq.isServed();   

  return pl.getChainComm().allReduceLand(hasReceived);
}


bool PendingSwap::isFinished()  
{
  auto result = _recvReq.isServed ();
  for(auto &elem : _sentReqs)
    result &= elem.isServed();
  return result; 
}


#else 

#endif
