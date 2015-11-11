#include "CredibleSet.hpp"
#include <cstring>
#include <unistd.h>
#include <cassert>
#include <sstream>
#include <iosfwd>

int NUM_BRANCHES; 

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS


static void printUsage(std::ostream &out )
{
  out << "credibleSet is a utility to extract a credible set of trees from a tree set.\n\n" ; 
  out << "USAGE: ./credibleSet -n id -f treeFile[..] [ -c credibleInterval]\n\n"; 
  out << "           -n id         a run id for the output file\n"; 
  out << "           -c int        credible set percentile (in (0,100], default:50)\n"; 
  out << "           -f file[..]   one or many topology input file(s)\n"; 
  out << std::endl; 
}


static auto  processCommandLine(int argc, char **argv)
  -> std::tuple<std::string,std::vector<std::string>,nat>
{
  int c = 0; 
  auto id = std::string{}; 
  nat credSet = 50; 
  auto files = std::vector<std::string>{}; 

  while (( c = getopt(argc, argv, "n:c:f:h"))   != EOF) 
    {
      switch(c)
	{
	case 'h': 
	  {
	    printUsage(std::cerr); 
	    exit(-1); 
	  }
	  break; 
	case 'n':
	  {
	    id = std::string{optarg, optarg + strlen(optarg)};
	  }
	  break; 
	case 'c' : 
	  {
	    auto &&iss = std::istringstream{optarg}; 
	    iss >> credSet ; 
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
	default: 
	  {
	    std::cerr << "Error: unknown option >" << char (c) << "<" << std::endl; 
	    exit(-1); 
	  }
	}
    }

  
  if(not ( 0 < credSet && credSet <= 100 ))
    {
      std::cerr << "Error: the percentile (passed via -c) must be between 0 (excluded) and 100 (included)." << std::endl; 
      exit(-1); 
    }
  
  return std::make_tuple(id,files,credSet);
}


int main(int argc, char **argv)
{
  NUM_BRANCHES = 1;   // BAD

  if(argc < 2 )
    {
      printUsage(std::cerr);
      exit(-1);       
    }

  auto id = std::string(argv[1]); 
  auto ciNumber =  0; 
  auto files = std::vector<std::string>{}; 

  std::tie(id,files,ciNumber) = processCommandLine(argc, argv);

  double ci = ciNumber; 
  ci /= 100.; 

  auto cs = CredibleSet(files); 

  auto &&ss = std::stringstream{}; 
  ss << PROGRAM_NAME << "_credibleSet." << id  ; 

  if(std::ifstream(ss.str()))
    {
      std::cerr << "The file " << ss.str() << " already exists (possibly left over from a previous run). Please\n"
		<< "choose a different run-id." << std::endl; 
      exit(-1); 
    }

  cs.printCredibleSet(ss.str(), ci); 

  return 0; 
}
