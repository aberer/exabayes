#ifndef _ABSTRACT_PENDING_SWAP_HPP
#define _ABSTRACT_PENDING_SWAP_HPP

#include "comm/SwapElem.hpp"
#include "comm/CommRequest.hpp"
#include "comm/PendingSwap.hpp"

class AbstractPendingSwap
{
public: 
  AbstractPendingSwap(SwapElem swap )
    : _swap(swap)
  {
  } 
  
  SwapElem getSwap () const {return _swap; } 

  virtual std::vector<char> getRemoteData() const = 0    ; 
  virtual bool isFinished() = 0 ;   
  virtual bool allHaveReceived(ParallelSetup& pl)  = 0 ; 
  virtual void initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid)  = 0; 

protected: 

  uint64_t cantorPair(uint64_t a, uint64_t b ) const 
  {
    return (a + b ) * (a + b + 1 ) / 2 + b ; 
  }

protected: 

  SwapElem _swap; 
  

}; 

#endif
