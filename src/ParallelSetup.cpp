#include "ParallelSetup.hpp"
#include "GlobalVariables.hpp"

#include <iostream> 


ParallelSetup::ParallelSetup(int argc, char **argv)
  : myRunBatch(0)
  , runsParallel(1)
  , globalRank(0)
  , globalSize(1)
{
#if HAVE_PLL == 0 
  MPI_Init(&argc, &argv);
#endif
}

#if HAVE_PLL == 0 
void ParallelSetup::initializeExaml(const CommandLine &cl)
{  
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  MPI_Comm_size(MPI_COMM_WORLD, &globalSize);  

  runsParallel = cl.getNumRunParallel();

  if(globalRank == 0 )
    std::cout << std::endl << std::endl << "This is " << PROGRAM_NAME << " process number: " << globalRank << " / " << globalSize << std::endl; 
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(globalSize <  runsParallel)
    {
      if(globalRank == 0)
	{
	  std::cout << "You requested to run "  << runsParallel << " in parallel, however there are only " 
		    << globalSize << " processes (we need at least 1 process per run, see command line option -R )"  << std::endl; 
	}
      MPI_Abort(MPI_COMM_WORLD, 1);
    }  

  if(globalRank == 0)
    std::cout << std::endl << runsParallel <<  " runs will be run in parallel." << std::endl; 

  // comm is the communicator used by the legacy axml-stuff in order to compute the likelihood.  
  int processesPerBatch = globalSize / runsParallel; 
  int myColor = globalRank / processesPerBatch; 
  int newRank = globalRank  % processesPerBatch; 
  
  myRunBatch = myColor; 
  // cout << "my color is " << myRunBatch <<  <endl; 


  MPI_Comm_split(MPI_COMM_WORLD, myColor, newRank, &comm); 
  
  MPI_Comm_rank(comm, &processID); 
  MPI_Comm_size(comm, &processes); 
  
  std::cout << std::endl <<  "\t\tprocess with global id "<< globalRank  << " works on batch "  << myRunBatch << " and has new rank " << processID << std::endl; 

}
#endif


void ParallelSetup::finalize()
{
#if HAVE_PLL == 0
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
}
