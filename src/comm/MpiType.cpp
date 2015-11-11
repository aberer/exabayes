
#include "common.h"

#if HAVE_PLL == 0 

#include "MpiType.hpp"


MPI_Datatype mpiType<int>::value = MPI_INT; 
MPI_Datatype mpiType<double>::value = MPI_DOUBLE; 
MPI_Datatype mpiType<char>::value = MPI_CHAR; 
MPI_Datatype mpiType<unsigned char>::value = MPI_UNSIGNED_CHAR; 

#else 
static void dummy()
{
}
#endif
