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


  int getMyRunBatch() const {return myRunBatch ; } 
  int getRunsParallel() const {return runsParallel ; }

private: 
  int myRunBatch; 
  int runsParallel;
  int globalRank; 		// 
  int globalSize; 
  
}; 



#endif
