#include "CommRequest.hpp"

#if HAVE_PLL == 0 		

#include <cassert>
#include <iostream>

#include "Communicator.hpp"
#include "extensions.hpp"


CommRequest::CommRequest( )
  : _req{nullptr}
{
}


CommRequest::CommRequest(CommRequest &&rhs)
  : _req(std::move(rhs._req))
  , _array(std::move(rhs._array))
{
}


void CommRequest::initialize( bool sending, int srcDest, int tag, std::vector<char> data,  Communicator& comm)
{
  assert(_req.get() == nullptr); 

  _req = make_unique<MPI_Request>();
  _array = data; 

  if(sending)
    MPI_Isend( (void*) _array.data(), _array.size(), MPI_CHAR, srcDest, tag, comm.getHandle(), _req.get());
  else 
    MPI_Irecv( (void*) _array.data(), _array.size(), MPI_CHAR, srcDest, tag, comm.getHandle(), _req.get());
}

CommRequest& CommRequest::operator=(CommRequest &&rhs)
{
  swap(*this, rhs); 
  return *this; 
}



void swap(CommRequest &lhs, CommRequest& rhs) 
{
  std::swap(lhs._req, rhs._req);
  std::swap(lhs._array, rhs._array);
}


bool CommRequest::isServed()
{
  assert(_req != nullptr ); 

  if( _req != nullptr )
    {
      int flag = 0; 
      MPI_Test(_req.get(), &flag, MPI_STATUS_IGNORE);
      return flag != 0 ; 
    }
  else 
    {
      std::cout << "DANGER request was empty "  << std::endl; 
      return true; 
    }
}

#else 

static void dummyFun()
{
  // avoids warnings 
}

#endif 
