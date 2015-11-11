#include <sstream>
#include <cstring>
#include <unistd.h>
#include <iosfwd>

int NUM_BRANCHES; 	

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "ConsensusTree.hpp"


auto processCommandLine(int argc, char **argv)
  -> std::tuple<std::string,std::vector<std::string>,nat, nat,bool>
{
  auto id = std::string{}; 
  auto thresh = double(50); 
  bool isRefined = true; 
  auto files = std::vector<std::string>{}; 
  nat burnin = 0; 

  int c = 0; 
  while( ( c = getopt(argc, argv, "n:t:f:") ) != EOF )
    {
      switch(c)
	{
	case 'n': 
	  {
	    id = std::string{optarg, strlen(optarg)}; 
	  }
	  break; 
	case 't' : 
	  {
	    if(std::string{"MRE"}.compare(optarg) == 0 
	       || std::string{"mre"}.compare(optarg) == 0)
	      {
		thresh = 50; 
		isRefined = true; 
	      }
	    else 
	      {
		auto &&iss = std::istringstream{optarg}; 
		iss >> thresh	; 
		isRefined = false; 
	      }
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
	    auto &&iss = std::istringstream {optarg}; 
	    iss >> burnin; 
	  }
	  break; 
	default: 
	  {
	    std::cerr << "Unrecognized option >" <<  optarg << "<"  << std::endl; 
	    exit(-1); 
	  }
	}
    }

  if(thresh < 50 || thresh > 100 )
    {
      std::cerr << "error: correct values for -t in [50,100] or MRE." << std::endl; 
      exit(-1); 
    }

  if(files.size() == 0 )
    {
      std::cerr << "Please specfiy tree input files via -f" << std::endl; 
      exit(-1); 
    }
  if(id.compare("") == 0)
    {
      std::cerr << "Please specify a runid for the output via -n. " << std::endl; 
      exit(-1); 
    }
 
  return std::make_tuple(id, files, burnin, thresh, isRefined);
}



static void printUsage(std::ostream &out)
{
  out << "consense computes various flavours of consensus trees from sets of trees.\n\n"; 
  out << "Usage: ./consense -n id  -f file[..] [-t threshold] [-b burnin] \n\n" ; 
  out << "        -n runid         an id for the output file\n" ; 
  out << "        -t thresh        a threshold for the consenus tree. Valid values:\n"
      << "                         values between 50 (majority rule) and 100 (strict) or MRE\n" 
      << "                         (the greedily refined MR consensus).  Default: MRE\n" ; 
  out << "        -b burnin        number of trees to discard for each file (from start). Default: off\n"; 
  out << "        -f file[..]      one or more exabayes topology files\n\n" ;
}


int main(int argc, char **argv)
{
  NUM_BRANCHES = 1;   // BAD

  if(argc < 2 ) 
    {
      printUsage(std::cout);
      exit(-1); 
    }

  bool isMre = false;   
  auto threshold = double{0.}; 
 
  std::string id= {}; 
  nat burnin = 0; 
  auto files = std::vector<std::string>{}; 

  std::tie(id, files, burnin, threshold, isMre) = processCommandLine(argc, argv);

  for(auto &file : files)
    {
      if(not std::ifstream(file))
	{
	  std::cerr << "Error: could not open file >" << file << "<" << std::endl; 
	  exit(-1); 
	}    
    }

  assert(threshold > 1); 
  threshold /= 100.; 

  auto ct = ConsensusTree(files); 
  auto result = ct.getConsensusTreeString(threshold,isMre);

  auto&& ss = std::stringstream{}; 
  ss << PROGRAM_NAME << "_consensusTree." << id; 

  if(std::ifstream(ss.str()))
    {
      std::cerr << std::endl << "File " << ss.str() << " already exists (probably \n"
		<< "left over from a previous run). Please choose a new run-id or remove\n"
		<< "previous output files." << std::endl; 
      exit(-1); 
    }

  std::ofstream outfile(ss.str()); 
  outfile << result << std::endl; 

  std::cout << "Printed consensus tree to " << ss.str() << std::endl; 

  return 0; 
}
