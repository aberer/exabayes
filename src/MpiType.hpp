#ifndef _MY_123_MPI_TYPE_HPP
#define _MY_123_MPI_TYPE_HPP

#include <mpi.h>

template<typename T>  
struct mpiType 
{
  static MPI_Datatype value ;   // = MPI_CHAR
}; 


template<> 
struct mpiType<int>  
{
  static MPI_Datatype value; 
}; 



template<>
struct mpiType<double>
{
  static MPI_Datatype value;  //  = MPI_DOUBLE
}; 


template<>
struct mpiType<unsigned int>
{
  static MPI_Datatype value;  //  = MPI_UNSIGNED
}; 

template<>
struct mpiType<char>
{
  static MPI_Datatype value;  //  = MPI_CHAR
}; 



template<>
struct mpiType<unsigned char>
{
  static MPI_Datatype value;  //  = MPI_UNSIGNED_CHAR
}; 


#endif






