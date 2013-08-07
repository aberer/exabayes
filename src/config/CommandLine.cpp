#include <unistd.h>
#include <iostream>
#include <cstring>
#include "CommandLine.hpp"
#include "GlobalVariables.hpp"
#include "ParallelSetup.hpp"
#include "file/OutputFile.hpp"
#include "MemoryMode.hpp"


CommandLine::CommandLine(int argc, char **argv)
  : configFileName("")
  , alnFileName("")
  , runid("")
  , treeFile("")
  , workDir("")
  , runNumParallel(1)
  , chainNumParallel(1)
  , checkpointId("")
  , memoryMode(MemoryMode::RESTORING)
  , perPartitionDataDistribution(false)
{
  seed.v[0] = 0; 
  seed.v[1] = 0; 
  
  parse(argc, argv) ;
}


void CommandLine::printVersion(bool toInfofile )
{   
  (toInfofile ? tout : std::cout )  << "This is " << PROGRAM_NAME << ", version " << PACKAGE_VERSION << std::endl << std::endl
				    << "For bugs reports and feature inquiries, please send an email to " << PACKAGE_BUGREPORT << std::endl; 
}



void CommandLine::printHelp()
{
  printVersion(false); 

  std::cout << std::endl << "./exabayes -f binFile -s seed -n id [options..] "
	    << std::endl; 

  std::cout << "\n\n"
	    << "Mandatory Arguments: \n"
	    << "    -f binFile       a binary alignment file that has been created by the appropriate parser before (see manual for help)\n"
	    << "    -s seed          a master seed for the MCMC\n"
	    << "    -n ruid          a run id\n" 
	    << std::endl;     

  std::cout << "\n" 
	    <<  "Options:\n" 
	    << "    -v               print version and quit\n"
	    << "    -h               print this help\n" 
	    << "    -c confFile      a file configuring your " << PROGRAM_NAME << " run. For a template see the examples/ folder\n"
	    << "    -w dir           specify a working directory for output files\n"
	    << "    -r id            restart from checkpoint. Just specify the id of the previous run here. \n"
	    << "                       Make sure, ExaBayes can access all files from this previous run.\n"
	    << "    -R num           the number of runs (i.e., independent chains) to be executed in parallel\n"
	    << "    -C num           number of chains (i.e., coupled chains) to be executed in parallel\n"
	    << "    -Q               per-partition data distribution (use this only with many partitions, check manual\n"
	    << "                       for detailed explanation)\n"
	    << "    -M mode          specifies the memory versus runtime trade (NOT IMPLEMENTED)\n"
	    << "            0          fastest\n" 
	    << "            1          standard\n"
	    << "            2          TODO\n"
	    << std::endl; 

  ParallelSetup::genericExit(-1); 
}


void CommandLine::assertFileExists(std::string filename)
{
  FILE *fh = myfopen(filename.c_str(), "r");

  if(fh == NULL )
    {
      fclose(fh); 
      std::cerr << "could not file file " << filename << ". Aborting." << std::endl; 
      ParallelSetup::genericExit(-1); 
    }
  fclose(fh); 
}


/** 
    @brief parses the command line 
 */ 
void CommandLine::parse(int argc, char *argv[])
{
  int c ; 

  // TODO threads/ processes? 
  
  while( (c = getopt(argc,argv, "c:f:vhn:w:s:t:R:r:M:C:Q")) != EOF)
    {
      try
	{	  
	  switch(c)
	    {
	    case 'c': 		// config file 	  
	      {
		configFileName = std::string(strdup(optarg)); 
		assertFileExists(configFileName);
	      }
	      break; 
	    case 'f': 		// aln file 
	      alnFileName = std::string(strdup(optarg)); 
	      assertFileExists(alnFileName); 
	      break; 
	    case 'v':  		// version 
	      printVersion(false );
	      ParallelSetup::genericExit(0); 
	      break; 
	    case 'h': 		// help 
	      printHelp();
	      break; 
	    case 'n': 		// runid 
	      runid = std::string(optarg); 	  
	      break; 
	    case 't': 		// trees -- have that in the config file? 
	      treeFile = std::string(strdup(optarg)); 
	      break; 
	    case 'w':		// working dir  
	      workDir = std::string(strdup(optarg)); 
	      break; 
	    case 's': 		// seed 
	      seed.v[0] = std::stoi(optarg);
	      break; 
	    case 'r': 
	      checkpointId = strdup(optarg);   
	      break; 
	    case 'M': 
	      memoryMode = MemoryMode(std::stoi(optarg)); 
	      break; 
	    case 'C': 
	      chainNumParallel = std::stoi(optarg); 
	      break; 
	    case 'R': 
	      runNumParallel = std::stoi(optarg);
	    case 'Q': 
	      perPartitionDataDistribution = true; 
	      break; 	  
	    default: 
	      {
		std::cerr << "Encountered unknown command line option " <<  c 
			  << "\n\nFor an overview of program options, please use -h" << std::endl ; 
		// TODO mpi-finalize stuff 
		abort();
	      }
	    }
	}
      catch(const std::invalid_argument& ia)
	{
	  std::cerr << "Invalid argument >" << optarg << "< to option >" << reinterpret_cast<char*>(&c) << "<" << std::endl; 
	  ParallelSetup::genericExit(-1);
	}
    }  
  
  if(runid.compare("") == 0 )
    {
      std::cerr << "please specify a runid with -n runid" << std::endl; 
      abort(); 
    }

  if(seed.v[0] == 0 && not checkpointId.compare("") == 0 )
    {
      std::cerr << "please specify a seed via -s seed (must NOT be 0)"   << std::endl; 
      abort(); 
    }


  if(seed.v[0] != 0 && checkpointId.compare("") != 0 )
    {
      std::cout << std::endl << "You provided a seed and run-id for a restart from a checkpoint.\n"
		<< "Please be aware that the seed will be ignored." << std::endl; 
    }

  
#if HAVE_PLL != 0
  if(runNumParallel > 1 || chainNumParallel > 1 )
    {
      std::cout << std::endl << "Your command line indicates that you intend to execute multiple runs\n"
		<< "or chains in parallel. This is the sequential version of" << PROGRAM_NAME << "\n"
		<< "and thus these command line flags will be ignored." << std::endl; 
    }

  if(perPartitionDataDistribution)
    {
      std::cout << std::endl << "Noticed your intent to enable per-partition data distribution (-Q). \n"
		<< "Since this is the sequential version of " << PROGRAM_NAME << "this option is ignored.\n"
		<< "For a detailed explanation of the -Q option, please consult the manual." << std::endl; 
    }
#endif



  if(workDir.compare("") != 0 && not OutputFile::directoryExists(workDir))
    {
      std::cout << std::endl << "Could not find the provided working directory >" << workDir << "<" << std::endl; 
      ParallelSetup::genericExit(-1);
    }

  if(alnFileName.compare("") == 0 )
    {
      std::cerr << "please specify a binary alignment file via -f file" <<  std::endl 
		<< "You have to transform your NEWICK-style alignment into a binary file using the appropriate parser (see manual)." << std::endl; 
      abort();
    }

}


randCtr_t CommandLine::getSeed() const
{
  return seed ; 
} 
