#include "Communicator.hpp"
#include <algorithm>
#include <numeric>

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

template<typename T>
std::vector<T> Communicator::gather(std::vector<T> myData, nat root) const
{
  auto result = std::vector<T>{};  
  nat totalLength = getSize() * myData.size(); 
  result.resize(totalLength); 

  MPI_Gather( myData.data(), myData.size(), mpiType<T>::value, result.data(), myData.size(), mpiType<T>::value, root, _comm);

  return result; 
}

// list all instantiations 
template std::vector<char> Communicator::gather(std::vector<char> myData, nat root) const; 
template std::vector<int> Communicator::gather(std::vector<int> myData, nat root) const; 


// could also return vector<vector<T>>...however, functionality not needed right now 
template<typename T> std::vector<T> Communicator::gatherVariableLength(std::vector<T> myData, int root) const 
{
  // determine lengths 
  auto lengths = std::vector<nat>{}; 
  auto myLen = std::vector<int>{  int(myData.size()) }; 
  auto allLengths = gather<int>( myLen, root );
  
  // calculate the displacements 
  auto displ = std::vector<int>{}; 
  displ.push_back(0);
  for(nat i = 1; i < allLengths.size(); ++i)
    displ.push_back(displ.back() + allLengths.at(i-1));

  auto result = std::vector<T>{}; 
  result.resize(std::accumulate(begin(allLengths), end(allLengths), 0));

  MPI_Gatherv(myData.data(), myData.size(), mpiType<T>::value, result.data(), allLengths.data(), displ.data(), mpiType<T>::value, root, _comm);
  
  return result; 
} 

template std::vector<char> Communicator::gatherVariableLength(std::vector<char> myData, int root) const ; 


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


template<typename T> 
T Communicator::receive( int source, int tag ) 
{
  auto result = T{}; 
  MPI_Recv(&result,1, mpiType<T>::value, source, tag, _comm, MPI_STATUS_IGNORE);
  return result; 
}

template<typename T> 
void Communicator::send( T elem, int dest, int tag ) 
{
  MPI_Send(&elem,1, mpiType<T>::value, dest, tag, _comm);
}


template int Communicator::receive(int source , int tag); 
template void Communicator::send( int elem, int dest, int tag ) ; 



#else 
#endif


