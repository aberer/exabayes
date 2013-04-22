/** 
    @file exabayes.cpp
    
    @brief This file sets the flavour of exabayes that has been
    compiled (i.e., sequential, pll, examl, ...)
    
*/ 



// TODO re-activate all that initial bla bla when starting up
// exa-bayes (model / program info )


#ifdef HAVE_AVX
#define __AVX
#endif

#include "axml.h" 
#include "bayes.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "main-common.h"
#include "proposals.h"
#include "output.h"
#include "adapters.h"

#include "CommandLine.hpp"

#if HAVE_PLL != 0

#include "globalVariables.h"

int main (int argc, char *argv[])
{ 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  assert(0); 
#endif

  ignoreExceptionsDenormFloat(); 

  CommandLine cl(argc, argv); 
  makeFileNames(); 
  exa_main(cl.getAdef(), cl.getSeed()); 

  return 0;
}


#else 

extern int processID; 
extern int processes;

int main(int argc, char *argv[])
{ 
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  printf("\nThis is %s FINE-GRAIN MPI Process Number: %d\n", PROGRAM_NAME, processID);   
  MPI_Barrier(MPI_COMM_WORLD);

  ignoreExceptionsDenormFloat(); 

  CommandLine cl(argc, argv); 

  if(processID == 0)
    makeFileNames(); 
  exa_main(cl.getAdef(), cl.getSeed()); 

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}


#endif




