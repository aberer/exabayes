#ifndef _LOCAL_SWAP_HPP
#define _LOCAL_SWAP_HPP

#include "comm/AbstractPendingSwap.hpp" 

#include <iostream>


class LocalSwap :  public AbstractPendingSwap
{
public: 
  LocalSwap( SwapElem elem )
    : AbstractPendingSwap(elem)
    , _haveReceived(false)
  { }

  virtual std::vector<char> getRemoteData()  const ; 
  virtual bool isFinished()  ;   
  virtual bool allHaveReceived(ParallelSetup& pl)   ; 
  virtual void initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid)  ; 

private: 
  // slight hack, because i was too lazy to change the signatures...  
  ParallelSetup* _plPtr; 
  
  bool _haveReceived; 
  std::vector<char> _dataReceived; 
  nat _runid; 
}; 

#endif
