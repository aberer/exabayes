#include <gtest/gtest.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#ifdef _WITH_MPI
#include <mpi.h>
#endif

#include "TopLevelInvocation.cpp" // 


int main (int argc, char **argv)
{
  int result = 0; 

  ::testing::InitGoogleTest(&argc, argv); 
#ifdef _WITH_MPI  
  MPI_Init(&argc, &argv); 
  
  int myRank = 0;  
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if(myRank > 0 )
    {
      // hack to silent remaining processes  
      int bak, newOne; 
      fflush(stdout); 
      bak = dup(1); 
      newOne = open("/dev/null", O_WRONLY) ; 
      dup2(newOne,1); 
      close(newOne); 
    }

  result = RUN_ALL_TESTS(); 
  MPI_Finalize(); 
#else 
  result = RUN_ALL_TESTS(); 
#endif
  
  return result; 
}
