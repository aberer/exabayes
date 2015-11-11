
#include "common.h"

#if HAVE_PLL == 0 

#include "MpiType.hpp"

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

MPI_Datatype mpiType<int>::value = MPI_INT; 
MPI_Datatype mpiType<double>::value = MPI_DOUBLE; 
MPI_Datatype mpiType<char>::value = MPI_CHAR; 
MPI_Datatype mpiType<unsigned char>::value = MPI_UNSIGNED_CHAR; 
MPI_Datatype mpiType<unsigned int>::value = MPI_UNSIGNED; 

#else 
static void dummy()
{
}
#endif
