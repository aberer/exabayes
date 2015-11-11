#ifndef _COMMUNICATOR_HPP
#define _COMMUNICATOR_HPP


#include "common.h"

#if HAVE_PLL == 0

#include <vector>
#include <mpi.h>
#include "MpiType.hpp"

class Communicator
{
public: 
  Communicator() ; 
  Communicator(const Communicator &rhs) = delete ; 
  Communicator(Communicator&& rhs); 
  Communicator& operator=( Communicator rhs); 
  friend void swap(Communicator &lhs, Communicator& rhs); 
  ~Communicator(); 

  void waitAtBarrier(); 

  MPI_Comm* getPtr(); 
  MPI_Comm& getHandle(); 

  int getRank() const ; 
  int getSize() const ; 
  bool isValid() const; 

  bool allReduceLand( bool myValue ) const ; 
  bool broadcast(bool value, int root) const ; 
  
  Communicator split(int color, int rank); 
  
  template<typename T> std::vector<T> gather(std::vector<T> myData, nat root = 0) const;
  template<typename T> std::vector<T> gatherVariableLength(std::vector<T> myData, int root = 0) const ; 
  
  template<typename T> T receive( int source, int tag ); 
  template<typename T> void send( T elem, int dest, int tag ) ; 

  template<typename T> std::vector<T> allReduce( std::vector<T> myValues); 

private: 
  MPI_Comm _comm; 
}; 


#else
class Communicator
{
}; 
#endif

#endif


