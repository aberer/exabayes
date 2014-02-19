#ifndef _REMOTE_COMM_IMPL_SEQ_HPP
#define _REMOTE_COMM_IMPL_SEQ_HPP

#include <vector>

#include "common.h"

class CommRequest; 

#include "comm/RemoteComm.hpp"

class RemoteComm::Impl
{
  typedef Impl SELF; 

public: 
  Impl(){} ; 
  Impl(const Impl& rhs)  {}; 
  Impl& operator=( Impl rhs); 
  friend void swap(Impl &lhs, Impl& rhs){}; 
  ~Impl(){}; 

  void createSendRequest(std::vector<char> array, int dest, int tag, CommRequest& req);
  void createRecvRequest(int src, int tag, nat length, CommRequest& req);  

  void waitAtBarrier() const {} 

#include "comm/CommCore.hpp"

  nat getNumberOfPhysicalNodes();
  static uint64_t _maxTagValue; 
  
}; 

#include "comm/dummy/RemoteCommImpl.tpp"


#endif
