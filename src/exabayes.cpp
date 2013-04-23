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
#include "chain.h"
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

  initParamStruct *initParams = (initParamStruct*)exa_calloc(1,sizeof(initParamStruct));   
  parseConfigWithNcl(configFileName, &initParams);
  setupGlobals(initParams); 

  exa_main(cl.getAdef(), cl.getSeed(), initParams); 

  return 0;
}


#else 

extern int processID; 
extern int processes;
extern MPI_Comm comm; 



int main(int argc, char *argv[])
{ 
  int globalCommSize = 0, 
    globalId = 0; 
  
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &globalId);
  MPI_Comm_size(MPI_COMM_WORLD, &globalCommSize);

  printf("\nThis is %s FINE-GRAIN MPI Process Number: %d / %d\n", PROGRAM_NAME, globalId, globalCommSize);   
  MPI_Barrier(MPI_COMM_WORLD);

  ignoreExceptionsDenormFloat(); 

  CommandLine cl(argc, argv); 

  initParamStruct *initParams = (initParamStruct*)exa_calloc(1,sizeof(initParamStruct));   
  parseConfigWithNcl(configFileName, &initParams);
  setupGlobals(initParams); 

  // comm is the communicator used by the legacy axml-stuff in order to compute the likelihood.  
  int processesPerBatch = globalCommSize / initParams->numRunParallel; 
  int myColor = globalId / processesPerBatch; 
  int newRank = globalId  % processesPerBatch; 

  MPI_Comm_split(MPI_COMM_WORLD, myColor, newRank, &comm); 
  
  if(processID == 0)
    makeFileNames(); 
  exa_main(cl.getAdef(), cl.getSeed(), initParams); 

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}


#endif




