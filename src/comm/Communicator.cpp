#include "Communicator.hpp"
#include <algorithm>

#if HAVE_PLL == 0

#include "CommRequest.hpp"

Communicator::Communicator() 
{
  _comm = MPI_COMM_WORLD; 
}


Communicator::Communicator(Communicator &&rhs)
  : _comm(std::move(rhs._comm))
{
}


Communicator& Communicator::operator=( Communicator rhs)
{
  swap(*this, rhs);
  return *this; 
}
  
void swap(Communicator &lhs, Communicator& rhs)
{
  std::swap(lhs._comm, rhs._comm);
}

MPI_Comm* Communicator::getPtr()
{
  return &_comm; 
}


MPI_Comm& Communicator::getHandle()
{
  return _comm; 
}

  
Communicator::~Communicator()
{
  if( _comm != MPI_COMM_WORLD )
    MPI_Comm_free(&_comm);
}


int Communicator::getRank() const 
{
  int result = 0; 
  MPI_Comm_rank(_comm, &result);
  return result; 
}


int Communicator::getSize() const 
{
  int result = 0; 
  MPI_Comm_size(_comm, &result); 
  return result; 
}

  
bool Communicator::isValid() const
{
  return _comm != MPI_COMM_NULL; 
}


void Communicator::waitAtBarrier()
{
  MPI_Barrier(_comm);
} 


Communicator Communicator::split(int color, int rank)
{
  auto result = Communicator{}; 
  MPI_Comm_split(_comm, color, rank, result.getPtr());
  return result; 
}


// a simple gather
std::vector<char> Communicator::gather(std::vector<char> myData) const
{
  auto result = std::vector<char>{};  
  nat totalLength = getSize() * myData.size(); 
  result.resize(totalLength); 

  MPI_Gather( myData.data(), myData.size(), MPI_CHAR, result.data(), myData.size(), MPI_CHAR, 0, _comm);

  return result; 
}


bool Communicator::broadcast(bool value, int root) const 
{
  char myVal = value ? 1 : 0 ; 
  MPI_Bcast(&myVal, 1, MPI_CHAR, root, _comm);
  
  return myVal == 1 ?true : false; 
}

// NOT used at the moment, would be nice, if we did 
template<typename T>
std::vector<T> Communicator::allReduce( std::vector<T> myValues)
{
  MPI_Allreduce(MPI_IN_PLACE, myValues.data(), myValues.size(), mpiType<T>::value, MPI_SUM, _comm); 
  return myValues; 
}



// implement with templates once we have more allreduces   
bool Communicator::allReduceLand( bool myValue ) const 
{
  char val =  myValue ? 1 : 0 ; 
  
  MPI_Allreduce(MPI_IN_PLACE, &val, 1, MPI_CHAR, MPI_LAND, _comm);

  return val == 1 ? true : false; 
}



template<> int Communicator::receive(int source , int tag); 
template<> void Communicator::send( int elem, int dest, int tag ) ; 



#else 
#endif


