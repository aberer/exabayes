#include <iosfwd>
#include <cstring>
#include <unistd.h>

int NUM_BRANCHES ; 

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "BipartitionExtractor.hpp" 

static void printUsage(std::ostream &out)
{
  out << "extractBips is a utility for extracting bipartitions, branch lengths\n"
      << "associated with bipartitions and statistics about these branch lengths\n"
      << "from tree sets.\n\n"; 
  out << "USAGE: ./extractBips -n  id -f file[..] [-b burnin]\n\n" ; 
  out << "Options:\n"; 
  out << "         -n id            an id for the output file\n"  ; 
  out << "         -f file[..]      one or more topology files\n"; 
  out << "         -b burnin        a number of trees to be discarded (beginning from\n"
      << "                          the start of the file). Default: 0\n\n"; 
}


static std::tuple<std::string, std::vector<std::string>,nat> processCommandLine(int argc, char **argv)
{
  nat burnin = 0; 
  auto files = std::vector<std::string> {}; 
  auto id = std::string{}; 

  int c = 0; 
  while(   (c = getopt(argc, argv, "n:f:b:h") ) != EOF )
    {
      switch(c)
	{
	  switch(c)
	    {
	    case 'h': 
	      {
		printUsage(std::cout); 
		exit(-1); 
	      }
	      break; 
	    case 'n': 
	      {
		id = std::string{optarg, optarg + strlen(optarg)}; 
	      }
	      break; 
	    case 'f': 
	      {
		nat index  = optind -1; 
		while(index < argc)
		  {
		    auto next = std::string{argv[index]}; 
		    ++index; 
		    if(next[0] != '-') 
		      files.push_back(next);
		    else 
		      {
			optind =  index -1 ; 
			break; 
		      }
		  }
	      }
	      break; 
	    case 'b': 
	      {
		auto &&iss = std::istringstream{optarg}; 
		iss >> burnin; 
	      }
	      break; 
	    default : 
	      {
		std::cerr << "Error: unknown option >" << char(c) << "<. " << std::endl; 
		exit(-1); 
	      }
	    }
	}
    }
  
  return std::make_tuple(id, files, burnin);
}

int main(int argc, char** argv)
{
  NUM_BRANCHES = 1;   // BAD

  if(argc < 2) 
    {
      printUsage(std::cout);
      exit(-1); 
    }

  auto files = std::vector<std::string>();
  nat burnin = 0; 
  auto id = std::string{};

  std::tie(id,files,burnin) = processCommandLine(argc, argv); 

  for(auto &file : files)
    {
      if(std::ifstream(file))
	{
	  std::cerr << "Error: could not open file >" <<  file  << "<" << std::endl; 
	  exit(-1); 
	}    
    }

  auto bipEx = BipartitionExtractor(files, false);
  bipEx.extractBips<true>(0);
  bipEx.printBipartitions(id);
  bipEx.printBipartitionStatistics(id); 
  bipEx.printFileNames(id);
  bipEx.printBranchLengths(id);

  return 0; 
}


