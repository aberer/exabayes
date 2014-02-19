#include "comm/LocalSwap.hpp"
#include "comm/ParallelSetup.hpp"


std::vector<char> LocalSwap::getRemoteData() const 
{
  assert(_haveReceived); 
  return _dataReceived; 
}

 
bool LocalSwap:: isFinished()  
{
  // a normal pending swap has to wait here, until the remote
  // processes have consumed his array. this does not apply here
  return true; 
}
  
bool LocalSwap::allHaveReceived(ParallelSetup& pl) 
{
  auto myId = _swap.getMyId(pl,_runid); 
  auto remoteId = _swap.getRemoteId(pl,_runid);

  auto tag = cantorPair(cantorPair(remoteId, myId), _runid); // NOTE: inverse of initialize 
  assert(tag < pl.getMaxTag()); 

  auto &localComm = pl.getRunComm().getLocalComm(); 
  
  if(not _haveReceived)
    {
      auto wasThere = false; 
      auto data = std::vector<char>();
      std::tie(wasThere, data) = localComm.readAsyncMessage<char>( tag ); 
      if(wasThere)
	{
	  _dataReceived = data; 
	  _haveReceived = true; 
	}
    }
  
  auto tmp = std::vector<int>{  _haveReceived ? 1 : 0   };
  tmp = localComm.allReduce(tmp); 
  
  return tmp[0] == localComm.size(); 
}

 
void LocalSwap::initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid)  
{
  _plPtr = &pl; 
  _runid = runid; 
  // post a message, that should be read by all threads belonging to
  // the other chain

  auto myId = _swap.getMyId(pl,runid); 
  auto remoteId = _swap.getRemoteId(pl, runid); 
  auto tag = cantorPair(cantorPair(myId, remoteId),runid); 
  assert(tag < pl.getMaxTag()); 

  if(pl.isChainLeader())
    pl.getRunComm().getLocalComm().postAsyncMessage( myChainSer, pl.getChainComm().size() , tag );
} 
