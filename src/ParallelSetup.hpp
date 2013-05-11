#ifndef _PARALLEL_SETUP_H
#define _PARALLEL_SETUP_H

#include "CommandLine.hpp"


class ParallelSetup
{
public: 
  ParallelSetup(int argc, char **argv);
  void finalize(); 

#if HAVE_PLL == 0 
  void initializeExaml(const CommandLine &cl); 
#endif

private: 
  
public: 
  int myRunBatch; 
  int runsParallel;
  int globalRank; 
  int globalSize; 
  
}; 



#endif
