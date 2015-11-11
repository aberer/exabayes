#include "parser/PhylipParser.hpp"
#include "ParallelSetup.hpp"

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <unistd.h>
#include <cstring>

int NUM_BRANCHES ; 

void helpMessage()
{
  std::cout << "\nparser produces a binary output file, that can be fed into\n"
	    << "ExaBayes/Yggdrasil. This is recommendable for large runs with hundreds\n"
	    << "of processes.\n\n" ; 

  std::cout << "./parser -s alnFile -q modelFile -n outputFile\n"
	    << "./parser -s alnFile -m {DNA|PROT} -n outputFile\n\n\n" 
	    << "     -s alnFile             a phylip-style alignment file\n\n"
	    << "     -q modelFile           a RAxML-style model file\n\n"
	    << "     -m dataType            specifies a datatype (either DNA or PROT) for a\n"
	    << "                              single partition. Not needed, when a\n"
	    << "                              model file is given.\n\n"
	    << "     -n outputfile          name of the output file\n\n"
	    << std::endl; 
}

static void myExit(int code)
{
  exit(code); 
}


int main(int argc, char **argv)
{
  if(argc < 2)
    {
      helpMessage(); 
      myExit(-1); 
    }
  
  auto alignmentFile = std::string(""); 
  auto modelFile = std::string(""); 
  auto singlePartitionModel = std::string(""); 
  auto outputFileName = std::string(""); 
  
  int c = 0 ; 
  while(  ( c  = getopt( argc, argv, "s:q:m:n:") )  != EOF ) 
    {
      try
	{
	  switch(c)
	    {
	    case 's': 
	      alignmentFile = std::string(optarg); 
	      break; 
	    case 'q': 
	      modelFile = std::string(optarg); 
	      break; 
	    case 'm':
	      singlePartitionModel = std::string(optarg); 
	      break; 
	    case 'n': 
	      outputFileName = std::string(optarg); 
	      break; 
	    default:
	      {
		std::cerr << "Encountered unknown command line option " <<  c 
			  << "\n\nFor an overview of program options, please use -h" << std::endl ; 
		myExit(-1); 
	      }
	    }
	}
      catch(const std::invalid_argument& ia )
	{
	  std::cerr << "Invalid argument >" << optarg << "< to option >" << reinterpret_cast<char*>(&c) << "<" << std::endl; 
	  myExit(-1); 
	}
    }

  // some validation
  if(outputFileName.compare("") == 0 )
    {
      std::cout << "Please specify an output file name via -n." << std::endl; 
      myExit(-1); 
    }

  if(outputFileName.compare("") != 0 && std::ifstream(outputFileName + ".binary"))
    {
      std::cout << "Error: output file "  <<  outputFileName  << ".binary already exists." << std::endl; 
      myExit(-1); 
    }

  if(modelFile.compare("") == 0 && singlePartitionModel.compare("") == 0)
    {
      std::cout << "Please specify either a model file via -q or a data type for a single\n"
		<< "partition (either DNA or PROT) via -m." << std::endl; 
      myExit(-1); 
    }

  if(modelFile.compare("") != 0 && not std::ifstream{modelFile})
    {
      std::cout << "Error: model file provided, but could not open model file >"  << modelFile << "<" << std::endl; 
      myExit(-1); 
    }

  
  if(not std::ifstream{alignmentFile})
    {
      std::cout << "Error: could not open alignment file >" << alignmentFile << "<" << std::endl; 
      myExit(-1); 
    }

  
  bool useSinglePartition =  modelFile.compare("") == 0; 

  auto parser = PhylipParser(alignmentFile, 
			     useSinglePartition ? singlePartitionModel : modelFile, 
			     not useSinglePartition);
  parser.parse(); 
  
  parser.writeToFile(outputFileName + ".binary"); 

  std::cout << "wrote binary alignment file to " << outputFileName << ".binary" << std::endl; 
  
  return 0 ;
}
