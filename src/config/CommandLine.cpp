#include <unistd.h>
#include <iostream>
#include <cstring>
#include "CommandLine.hpp"
#include "GlobalVariables.hpp"
#include "ParallelSetup.hpp"
#include "file/OutputFile.hpp"
#include "MemoryMode.hpp"
#include "FlagType.hpp"

CommandLine::CommandLine(int argc, char **argv)
  : configFileName("")
  , alnFileName("")
  , runid("")
  , treeFile("")
  , workDir("")
  , runNumParallel(1)
  , chainNumParallel(1)
  , checkpointId("")
  , memoryMode(MemoryMode::RESTORE_ALL)
  , perPartitionDataDistribution(false)
  , saveMemorySEV(false)
  , dryRun (false)
  , modelFile("")
  , singleModel("")
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

  std::cout << std::endl << "./exabayes   -f alnFile [ -q modelFile ] [ -m model ] [ -s seed | -r id ]  -n id [options..] "
	    << std::endl; 

  std::cout << "\n\n"
	    << "Mandatory Arguments: \n"
	    << "    -f alnFile       a alignment file (either binary and created by parser or plain-text phylip)\n"
	    << "    -s seed          a master seed for the MCMC\n"
	    << "    -n ruid          a run id\n" 
	    << "    -r id            restart from checkpoint. Just specify the id of the previous run (-n) here. \n"
	    << "                       Make sure, that all files created by the previous run are in the working directory.\n"
	    << "                       This option is not mandatory for the start-up, seed (via -s) will be ignored.\n"
	    << "    -q modelfile     a RAxML-style model file (see manual) for multi-partition alignments. Not needed \n"
	    << "                       with binary files.\n"
	    << "    -t treeFile      a file containing starting trees (in Newick format) for chains. If the file provides less\n"
	    << "                       starting trees than chains to be initialized, parsimony/random trees will be used for\n"
	    << "                       remaining chains. If a tree contains branch lengths, these branch lengths will be used\n"
	    << "                       as initial values.\n"
	    << "    -m model         indicates the type of data for a single partition non-binary alignment file\n" 
	    << "                       (valid values: DNA or PROT)\n"
	    << std::endl;     

  std::cout << "\n" 
	    <<  "Options:\n" 
	    << "    -v               print version and quit\n"
	    << "    -h               print this help\n" 
	    << "    -d               execute a dry-run. Procesess the input, but does not execute any sampling.\n"
	    << "    -c confFile      a file configuring your " << PROGRAM_NAME << " run. For a template see the examples/ folder\n"
	    << "    -w dir           specify a working directory for output files\n"
	    << "    -R num           the number of runs (i.e., independent chains) to be executed in parallel\n"
	    << "    -C num           number of chains (i.e., coupled chains) to be executed in parallel\n"
	    << "    -Q               per-partition data distribution (use this only with many partitions, check manual\n"
	    << "                       for detailed explanation)\n"
	    << "    -S               try to save memory using the SEV-technique for gap columns on large gappy alignments\n" 
	    << "                       Please refer to  http://www.biomedcentral.com/1471-2105/12/470\n" 
	    << "                       On very gappy alignments this option yields considerable runtime improvements. \n"
	    << "    -M mode          specifies the memory versus runtime trade-off (see manual for detailed discussion).\n"
	    << "                       <mode> is a value between 0 (fastest, highest memory consumption) and 3 (slowest,\n"
	    << "                       least memory consumption)\n"
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
  
  while( (c = getopt(argc,argv, "c:df:vhn:w:s:t:R:r:M:C:Qm:Sq:")) != EOF)
    {
      try
	{	  
	  switch(c)
	    {
	    case 'c': 		// config file 	  
	      {
		configFileName = std::string(optarg); 
		assertFileExists(configFileName);
	      }
	      break; 
	    case 'f': 		// aln file 
	      alnFileName = std::string(optarg); 
	      assertFileExists(alnFileName); 
	      break; 
	    case 'v':  		// version 
	      printVersion(false );
	      ParallelSetup::genericExit(0); 
	      break; 
	    case 'd': 
	      dryRun = true; 
	      break; 
	    case 'h': 		// help 
	      printHelp();
	      break; 
	    case 'n': 		// runid 
	      runid = std::string(optarg); 	  
	      break; 
	    case 't': 		// trees -- have that in the config file? 
	      treeFile = std::string(optarg); 
	      break; 
	    case 'w':		// working dir  
	      workDir = std::string(optarg); 
	      break; 
	    case 's': 		// seed 
	      seed.v[0] = std::stoi(optarg);
	      break; 
	    case 'r': 
	      checkpointId = std::string{optarg};   
	      break; 
	    case 'q':
	      modelFile = std::string{optarg}; 
	      break; 
	    case 'm': 
	      singleModel = std::string{optarg}; 
	      break; 
	    case 'S': 
	      saveMemorySEV = true; 
	      break; 
	    case 'M': 
	      memoryMode = MemoryMode(std::stoi(optarg)); 
	      break; 
	    case 'C': 
	      chainNumParallel = std::stoi(optarg); 
	      break; 
	    case 'R': 
	      runNumParallel = std::stoi(optarg);
	      break; 
	    case 'Q': 
	      perPartitionDataDistribution = true; 
	      break; 	  
	    default: 
	      {
		std::cerr << "Encountered unknown command line option " <<  c 
			  << "\n\nFor an overview of program options, please use -h" << std::endl ; 
		// TODO mpi-finalize stuff 
		exit(-1); 
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
      exit(-1); 
    }
  
  if(seed.v[0] != 0 && checkpointId.compare("") != 0 )
    {
      std::cout << std::endl << "You provided a seed and run-id for a restart from a checkpoint.\n"
		<< "Please be aware that the seed will be ignored." << std::endl; 
    }

  if(checkpointId.compare("") == 0 && seed.v[0] == 0 )
    {
      std::cerr << "please specify a seed via -s seed (must NOT be 0)"   << std::endl; 
      exit(-1); 
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
      std::cerr << "please specify an alignment file via -f file" <<  std::endl 
		<< "You have to transform your NEWICK-style alignment into a binary file using the appropriate parser (see manual)." << std::endl; 
      exit(-1); 
    }

  if(alnFileIsBinary())
    {
      if(singleModel.compare("") != 0 || modelFile.compare("") != 0 )
	{
	  std::cout << "Found binary alignment file. Additionally, you provided a model file\n"
	    "(-q) or specified a data type for a single partiton. This information\n"
	    "will be ignored.\n"; 
	  modelFile = ""; 
	  singleModel = ""; 
	}
    }
  else 
    {
      if(singleModel.compare("") == 0 && ( modelFile.compare("") == 0 || not std::ifstream(modelFile) )  )
	{
	  std::cout << "Found a phylip-style alignment file. However, you did not provide a\n"
	       << "model file (see -q, resp. it coul not be found) or a data type specification for a single\n"
	       << "partition (-m). Cannot proceed.\n" ; 
	  exit(-1); 
	}
    }
}

RunModes CommandLine::getTreeInitRunMode() const 
{
  auto runmodes = RunModes::NOTHING; 
  if(isPerPartitionDataDistribution())
    runmodes = runmodes | RunModes::PARTITION_DISTRIBUTION; 
  if(isSaveMemorySEV())
    runmodes = runmodes | RunModes::MEMORY_SEV; 
  return runmodes; 
} 



randCtr_t CommandLine::getSeed() const
{
  return seed ; 
} 


bool CommandLine::alnFileIsBinary() const
{
  auto&& in =std::ifstream (alnFileName, std::ios::binary); 

  assert(in) ; 

  auto fileId = std::string{"BINARY"} ; 
  char firstBytes[7]; 
  memset(firstBytes, '\0', 7 * sizeof(char)); 
  in.read(firstBytes, 6 * sizeof(char));
  auto readString = std::string(firstBytes); 

  bool result = readString.compare(fileId) == 0 ;
  // if(not result)
  //   {
  //     // std::cout << "wanted to read >" << fileId << "< instead got >" << readString << "<" << std::endl; 
      
  //   }
  // else 
  //   {
  //     // std::cout << "Determined alignment file to be binary" << std::endl; 
  //   }

  return result; 
}  


