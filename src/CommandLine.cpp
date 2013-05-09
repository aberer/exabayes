#include <unistd.h>
#include <iostream>
#include "CommandLine.hpp"
#include "GlobalVariables.hpp"

using namespace std; 


CommandLine::CommandLine(int argc, char *argv[])
{
  parse(argc, argv) ;
}


void CommandLine::printVersion()
{
#if HAVE_PLL == 0
  if(processID == 0)
#endif
    {
      cout << "This is " << PROGRAM_NAME << ", version " << PACKAGE_VERSION << endl << endl
	   << "For bugs reports and feature inquiries, please send an email to " << PACKAGE_BUGREPORT << endl; 
    }
}



void CommandLine::printHelp()
{
#if HAVE_PLL == 0
  if(processID != 0)
    exit(0); 
#endif


  printVersion(); 

  cout << "\n\n" 
       <<  "Options:\n" 
       << "\t-f binFile\t\ta binary alignment file that has been created by the appropriate parser before (see manual for help) [mandatory] \n"
       << "\t-c confFile\t\ta config file. For a template see the examples/ folder [mandatory]\n"
       << "\t-v\t\tprint version and exit\n"
       << "\t-h\t\tprint this help\n" 
       << "\t-w dir\t\tspecify a working directory for output files (NOT IMPLEMENTED)\n"
       << "\t-s seed\t\ta master seed for the MCMC [mandatory]\n"
       << "\t-n ruid\t\ta run id [mandatory]\n"
       << endl; 

  exit(0); 
}





/** 
    @brief parses the command line 
 */ 
void CommandLine::parse(int argc, char *argv[])
{
  int c ; 


  bool runidSet = false, 
  configFileSet = false, 
  alnFileSet = false, 
    seedSet = false; 

  // TODO threads/ processes? 
  // TODO implement working directory 
  
  while( (c = getopt(argc,argv, "c:f:vhn:w:s:t:")) != EOF)
    {
      switch(c)
	{
	case 'c': 		// config file 	  
	  {
	  strcpy(configFileName, optarg) ; 
	  if( NOT filexists(configFileName))
	    {
	      cout << "could not find file "  << configFileName << ". Aborting." << endl; 
	      exit(0); 
	    }
	  configFileSet = true; 
	  }
	  break; 
	case 'f': 		// aln file 
	  strcpy(byteFileName, optarg); 
	  if( NOT filexists(byteFileName))
	    {
	      cout << "could not find file "  << byteFileName << ". Aborting." << endl; 
	      exit(0); 
	    }
	  alnFileSet = true; 
	  break; 
	case 'v':  		// version 
	  printVersion();
	  exit(0); 
	  break; 
	case 'h': 		// help 
	  printHelp();
	  break; 
	case 'n': 		// runid 
	  strcpy(run_id,optarg); 
	  runidSet = true ; 
	  break; 
	case 't': 		// trees -- have that in the config file? 
	  strcpy(tree_file, optarg); 
	  if( NOT filexists(tree_file))
	    {
	      cout << "could not find file "  << tree_file << ". Aborting." << endl; 
	      exit(0); 
	    }
	  break; 
	case 'w':		// working dir TODO 
	  strcpy(workdir, optarg); 
	  break; 
	case 's': 		// seed 
	  seed = atoi(optarg);
	  seedSet = true; 
	  break; 
	default: 
	  printf("?? Encountered unknown command line argument 0%o ??\n\nFor an overview of program options, please use -h", c);
	}
    }  
  
  if(not runidSet)
    {
      cerr << "please specify a runid with -n runid" << endl; 
      abort(); 
    }

  if(not seedSet)
    {
      cerr << "please specify a seed via -s seed"   << endl; 
      abort(); 
    }

  if(not configFileSet)
    {
      cerr << "please specify a config file via -c configfile" << endl << "You find a template in the examples/ folder" << endl; 
      abort(); 
    }

  if(not alnFileSet)
    {
      cerr << "please specify a binary alignment file via -f file" <<  endl 
	   << "You have to transform your NEWICK-style alignment into a binary file using the appropriate parser (see manual)." << endl; 
      abort();
    }

}




