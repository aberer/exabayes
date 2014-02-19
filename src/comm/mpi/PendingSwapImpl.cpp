#include "PendingSwapImpl.hpp" 

#include "comm/ParallelSetup.hpp"

PendingSwap::Impl::Impl(SwapElem elem)
  : AbstractPendingSwap(elem)
  {
  }


std::vector<char> PendingSwap::Impl::getRemoteData()  const 
{
  assert(_recvReq.isServed()); 
  return _recvReq.getArray(); 
}


uint64_t PendingSwap::Impl::createTag( int rank ) const 
{
  auto chainA = uint64_t(_swap.getOne()); 
  auto chainB = uint64_t(_swap.getOther()); 
  auto result = cantorPair(uint64_t(rank), cantorPair(chainA,chainB)); 
  return result; 
}


void PendingSwap::Impl::initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid) 
{
  uint64_t tag = createTag( pl.getChainComm().getRank() ); 
  // std::cout << SyncOut() << SHOW(tag) << std::endl; 
  assert(tag < pl.getMaxTag()); 

  nat myId = _swap.getMyId(pl,runid); 
  nat remoteId = _swap.getRemoteId(pl,runid);

  nat ourRunBatch = runid % pl.getRunsParallel();
  nat myBatch = myId % pl.getChainsParallel();
  nat remoteBatch = remoteId % pl.getChainsParallel(); 
  
  nat numInMyBatch = pl.getMesh().getNumRanksInDim(ourRunBatch, myBatch); 
  nat numInRemoteBatch = pl.getMesh().getNumRanksInDim(ourRunBatch, remoteBatch); 

  auto myCoords = pl.getMyCoordinates(); 

  // determine from who i receive 
  // auto srcRankInRun =  ; 
  auto remRank =  pl.getRunComm().mapToRemoteRank( pl.getRankInRun( { {ourRunBatch, remoteBatch, myCoords[2] % numInRemoteBatch }  } )); 

  nat myRankInRun =  pl.getRankInRun( myCoords ); 

  assert( int(myRankInRun) == pl.getRunComm().getRank() ); 

  pl.getRunComm().createRecvRequest(remRank, tag, myChainSer.size(), _recvReq);

  // send stuff 
  for(nat i = 0 ; i < numInRemoteBatch; ++i)
    {
      if( i % numInMyBatch == myCoords[2])
	{
	  auto hisRemoteRank = pl.getRunComm().mapToRemoteRank( pl.getRankInRun( {{ ourRunBatch, remoteBatch, i } })); 
	  
	  _sentReqs.push_back(  CommRequest()  );
	  auto &&elem = _sentReqs.back();
	  pl.getRunComm().createSendRequest(myChainSer, hisRemoteRank, tag , elem); 
	}
    }
}


bool PendingSwap::Impl::allHaveReceived(ParallelSetup& pl)  
{
  bool hasReceived = _recvReq.isServed();   

  //  boolean operators would be nicer, but not every mpi
  // implementation offers these

  nat toReduce = hasReceived ? 1 : 0 ; 
  auto arr = std::vector<nat>{toReduce};
  arr = pl.getChainComm().allReduce(arr); 

  toReduce = arr.at(0); 
    
  // std::cout << SyncOut()  << pl.getGlobalComm().getRank() << " received by "<< toReduce << std::endl; 

  return toReduce == nat(pl.getChainComm().size()); 
}



bool PendingSwap::Impl::isFinished()  
{
  auto result = _recvReq.isServed ();
  for(auto &elem : _sentReqs)
    result &= elem.isServed();
  return result; 
}


