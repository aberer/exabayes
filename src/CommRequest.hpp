#ifndef _COMM_REQUEST_HPP
#define _COMM_REQUEST_HPP

#include "common.h"

#if HAVE_PLL == 0 

#include <vector>
#include <memory>
#include <mpi.h>
#include <cassert>

#include "Communicator.hpp"

class CommRequest
{
public: 
  CommRequest(); 
  CommRequest(CommRequest &&rhs); 
  CommRequest(const CommRequest& rhs) = delete; 
  CommRequest& operator=(CommRequest &&rhs)  ; 
  friend void swap(CommRequest &lhs, CommRequest& rhs) ; 

  void initialize( bool sending, int srcDest, int tag, std::vector<char> data, Communicator& comm);   
  std::vector<char> getArray() const {return _array; }
  bool isServed(); 

private: 
  std::unique_ptr<MPI_Request> _req; 
  std::vector<char> _array; 
}; 

#endif

#endif
