#include <unistd.h>
#include <iostream>

#include "CommandLine.hpp"
#include "GlobalVariables.hpp"

using namespace std; 


CommandLine::CommandLine(int argc, char *argv[])
  : seed(0)
  , configFileName("")
  , alnFileName("")
  , runid("")
  , treeFile("")
  , workDir("")
  , runNumParallel(1)
{
  parse(argc, argv) ;
}


void CommandLine::printVersion(bool toInfofile )
{   
  (toInfofile ? tout : cout )  << "This is " << PROGRAM_NAME << ", version " << PACKAGE_VERSION << endl << endl
    			       << "For bugs reports and feature inquiries, please send an email to " << PACKAGE_BUGREPORT << endl; 
}



void CommandLine::printHelp()
{
#if HAVE_PLL == 0
  if(processID != 0)
    exit(0); 
#endif

  printVersion(false); 

  cout << "\n\n" 
       <<  "Options:\n" 
       << "\t-f binFile\t\ta binary alignment file that has been created by the appropriate parser before (see manual for help) [mandatory] \n"
       << "\t-c confFile\t\ta config file. For a template see the examples/ folder [mandatory]\n"
       << "\t-v\t\tprint version and exit\n"
       << "\t-h\t\tprint this help\n" 
       << "\t-w dir\t\tspecify a working directory for output files (NOT IMPLEMENTED)\n"
       << "\t-s seed\t\ta master seed for the MCMC [mandatory]\n"
       << "\t-n ruid\t\ta run id [mandatory]\n"
       << "\t-R num\t\tthe number of runs (i.e., independent chains) to be executed in parallel"
       << endl; 

  exit(0); 
}


void CommandLine::assertFileExists(string filename)
{
  FILE *fh = myfopen(filename.c_str(), "r");

  if(fh == NULL )
    {
      fclose(fh); 
      cerr << "could not file file " << filename << ". Aborting." << endl; 
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
	    
	  }
	  break; 
	case 'c': 		// config file 	  
	  {
	    configFileName = string(optarg); 
	    assertFileExists(configFileName);
	  }
	  break; 
	case 'f': 		// aln file 
	  alnFileName = string(optarg); 
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
	  runid = string(optarg); 	  
	  break; 
	case 't': 		// trees -- have that in the config file? 
	  treeFile = string(optarg); 
	  break; 
	case 'w':		// working dir  
	  workDir = string(optarg); 
	  break; 
	case 's': 		// seed 
	  seed = atoi(optarg);
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
      cerr << "please specify a runid with -n runid" << endl; 
      abort(); 
    }

  if(seed == 0 )
    {
      cerr << "please specify a seed via -s seed (must NOT be 0)"   << endl; 
      abort(); 
    }

  if( configFileName.compare("") == 0)
    {
      cerr << "please specify a config file via -c configfile" << endl << "You find a template in the examples/ folder" << endl; 
      abort(); 
    }

  if(alnFileName.compare("") == 0 )
    {
      cerr << "please specify a binary alignment file via -f file" <<  endl 
	   << "You have to transform your NEWICK-style alignment into a binary file using the appropriate parser (see manual)." << endl; 
      abort();
    }

}




