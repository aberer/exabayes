#include <unistd.h>
#include <iostream>
#include <cstring>
#include "CommandLine.hpp"
#include "GlobalVariables.hpp"


CommandLine::CommandLine(int argc, char **argv)
  : configFileName("")
  , alnFileName("")
  , runid("")
  , treeFile("")
  , workDir("")
  , runNumParallel(1)
  , checkpointId("")
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
#if HAVE_PLL == 0
  if(processID != 0)
    exit(0); 
#endif

  printVersion(false); 

  std::cout << "\n\n" 
       <<  "Options:\n" 
       << "\t-f binFile\t\ta binary alignment file that has been created by the appropriate parser before (see manual for help) [mandatory] \n"
       << "\t-c confFile\t\ta config file. For a template see the examples/ folder [mandatory]\n"
       << "\t-v\t\tprint version and exit\n"
       << "\t-h\t\tprint this help\n" 
       << "\t-w dir\t\tspecify a working directory for output files (NOT IMPLEMENTED)\n"
       << "\t-s seed\t\ta master seed for the MCMC [mandatory]\n"
       << "\t-n ruid\t\ta run id [mandatory]\n"
       << "\t-R num\t\tthe number of runs (i.e., independent chains) to be executed in parallel"
	    << std::endl; 

  exit(0); 
}


void CommandLine::assertFileExists(std::string filename)
{
  FILE *fh = myfopen(filename.c_str(), "r");

  if(fh == NULL )
    {
      fclose(fh); 
      std::cerr << "could not file file " << filename << ". Aborting." << std::endl; 
      exit(0);
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
  // TODO implement working directory 
  
  while( (c = getopt(argc,argv, "c:f:vhn:w:s:t:R:")) != EOF)
    {
      switch(c)
	{
	case 'C': 		// chain-level parallelism 
	  {
	    std::cerr << "not in use" << std::endl; 
	    assert(0); 
	  }
	  break; 
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
	  exit(0); 
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
	  
	  
	  break; 
	case 'R': 
	  runNumParallel = atoi(optarg);
	  break; 	  
	default: 
	  {
	    printf("?? Encountered unknown command line argument 0%o ??\n\nFor an overview of program options, please use -h", c);
	    // TODO mpi-finalize stuff 
	    abort();
	  }
	}
    }  
  
  if(runid.compare("") == 0 )
    {
      std::cerr << "please specify a runid with -n runid" << std::endl; 
      abort(); 
    }

  if(seed.v[0] == 0 )
    {
      std::cerr << "please specify a seed via -s seed (must NOT be 0)"   << std::endl; 
      abort(); 
    }

  if( configFileName.compare("") == 0)
    {
      std::cerr << "please specify a config file via -c configfile" << std::endl << "You find a template in the examples/ folder" << std::endl; 
      abort(); 
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
