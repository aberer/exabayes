#include "ParallelSetup.hpp"
#include "GlobalVariables.hpp"

#include <iostream> 

using namespace std; 


ParallelSetup::ParallelSetup(int argc, char **argv)
{
#if HAVE_PLL == 0 
  MPI_Init(&argc, &argv);
#endif
}

#if HAVE_PLL == 0 
extern MPI_Comm comm; 		// the examl communicator 
extern int processID; 		// examl rank 
extern int processes; 		// examl comm size 

void ParallelSetup::initializeExaml(const CommandLine &cl)
{  
  runsParallel = cl.getNumRunParallel();

  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  MPI_Comm_size(MPI_COMM_WORLD, &globalSize);
  
  tout << endl << endl << "This is " << PROGRAM_NAME << " process number: " << globalRank << " / " << globalSize << endl; 
  MPI_Barrier(MPI_COMM_WORLD);

  
  if(globalSize <  runsParallel)
    {
      if(globalRank == 0)
	{
	  tout << "You requested to run "  << runsParallel << " in parallel, however there are only " 
	       << globalSize << " processes (we need at least 1 process per run, see command line option -R )"  << endl; 
	}
      MPI_Abort(MPI_COMM_WORLD, 1);
    }  

  // comm is the communicator used by the legacy axml-stuff in order to compute the likelihood.  
  int processesPerBatch = globalSize / runsParallel; 
  int myColor = globalRank / processesPerBatch; 
  int newRank = globalRank  % processesPerBatch; 
  
  MPI_Comm_split(MPI_COMM_WORLD, myColor, newRank, &comm); 

  MPI_Comm_rank(comm, &processID); 
  MPI_Comm_size(comm, &processes); 

  printf("\n\n process %d working on batch %d\n", processID, myColor); 
}
#endif


void ParallelSetup::finalize()
{
#if HAVE_PLL == 0
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
}
