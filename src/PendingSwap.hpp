#ifndef _UNFINISHED_SWAP
#define _UNFINISHED_SWAP

#include <cassert>
#include <vector>

#include "common.h"
#include "SwapElem.hpp"
#include "extensions.hpp"

/** 
    @file PendingSwap.hpp

    @brief represents a swap for which in the parallel case have not received the other chain yet

    This class has no meaning in non-MPI-versions. 
*/ 



#if HAVE_PLL == 0

#include <mpi.h>
#include "CommRequest.hpp"

class ParallelSetup; 

class PendingSwap
{
public: 
  PendingSwap( SwapElem swap); 
  PendingSwap(PendingSwap&& rhs); 
  PendingSwap& operator=( PendingSwap elem) = delete ; 
  PendingSwap( const PendingSwap& rhs) = delete; 

  friend void swap(PendingSwap& lhs, PendingSwap& rhs);

  bool isFinished();

  void initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid, nat numChains); 
  SwapElem getSwap() const {return _swap; }

  bool allHaveReceived(ParallelSetup& pl) ;

  std::string getRemoteData()  ; 
  
private: 			// METHODS 
  int createTag( nat numChains ) const ; 

private: 			// ATTRIBUTES
  std::vector<CommRequest> _sentReqs; 
  CommRequest _recvReq; 
  SwapElem _swap; 
}; 

#else 


#include <string>

class PendingSwap
{
public:   
  PendingSwap(SwapElem swap)
    : _swap{swap}
  {
  }

  bool isFinished()  
  {
    assert(0); 
    return true; 
  }


  std::string getRemoteData()  
  {
    assert(0);
    auto result = std::string{};
    return result; 
  } 



  bool allHaveReceived(ParallelSetup& pl)  
  {
    assert(0); 
    return true; 
  }

  void execute(){assert(0); }

  SwapElem getSwap() const { assert(0); return _swap; }

  void initialize(ParallelSetup& pl, std::vector<char> myChainSer, nat runid, nat numChains)
  {
    assert(0); 
  } 

private: 
  SwapElem _swap; 
}; 


#endif

#endif

