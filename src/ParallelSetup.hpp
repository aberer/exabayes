#ifndef _PARALLEL_SETUP_H
#define _PARALLEL_SETUP_H

#include "config/CommandLine.hpp"

#if HAVE_PLL == 0

extern MPI_Comm comm; 		// the examl communicator 
extern int processID; 		// examl rank 
extern int processes; 		// examl comm size 

#endif

#include <iostream>

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
  nat getProcessesPerBatch() const {return globalSize / runsParallel; }
  nat getMyIdOnRunLevel() const {return globalRank  % getProcessesPerBatch();  }
  bool isReportingProcess() const {  return globalRank == 0 ; }
  
private: 
  int myRunBatch; 
  int runsParallel;
  int globalRank; 		// 
  int globalSize; 
  
}; 



#endif
