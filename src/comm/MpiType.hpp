#ifndef _MY_123_MPI_TYPE_HPP
#define _MY_123_MPI_TYPE_HPP

#include <mpi.h>

template<typename T>  
struct mpiType 
{
  static MPI_Datatype value ;   // = MPI_CHAR
}; 


#endif






