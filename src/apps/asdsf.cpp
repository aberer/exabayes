#include <cassert>
#include <unistd.h>
#include <unordered_map>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "axml.h"

int NUM_BRANCHES; 

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "common.h"

#include "SplitFreqAssessor.hpp"


void myExit(int code)
{
  exit(code); 
}



void printUsage()
{
  std::cout
	    << "\nComputes the avgerage and maximum deviation of split frequencies\n"
	    << "of sets of trees.\n\n"
	    << "USAGE: ./asdsf [-m] [-b relBurnin | -r start ] [ -i ignoreFreq ]  -f file[..]\n\n"
	    << "      -m               if files contain a different number of trees, try \n"
	    << "                       to use as many trees as possible [default: use\n"
	    << "                       first n trees, where n is the minimum number\n"
	    << "                       available in all tree files]\n\n"  
	    << "      -i ignoreFreq    ignore splits with frequencies lower than ignoreFreq\n"
	    << "                       [Range: 0.0 < ignoreFreq < 1.0; default: 0.1]\n\n"
	    << "      -f file[..]      two or more topology files\n\n"
	    << "      -b relBurnin     discard first relBurnin percent of tree samples \n"
	    << "                       [ 0.0 <= relBurnin < 1.0; default 0.25]\n\n"
	    << "      -r num           constant burn-in: discard the first <num> samples\n"
	    << "\n\n"
    ; 

  myExit(0); 
}



int main(int argc, char** argv)
{

// #ifdef _USE_GOOGLE_PROFILER
//   auto myProfileFile = "profile.out"; 
//   remove(myProfileFile);
//   ProfilerStart(myProfileFile);
// #endif


  NUM_BRANCHES = 1; 

  bool constBurninWasSet = false; 
  bool relBurninWasSet = false; 

  nat constBurnin = 0; 
  double relBurnin = 0.25; 
  auto files = std::vector<std::string>{}; 
  double ignoreFreq = 0.1; 

  int c = 0;   
  bool takeAll = false; 
  while( ( c = getopt(argc, argv,"r:b:f:mi:h") ) != EOF )
    {
      switch(c)
	{
	case 'h': 
	  {
	    printUsage();
	  }
	  break; 
	case 'i': 
	  {
	    auto &&iss = std::stringstream{optarg}; 
	    iss >> ignoreFreq ; 
	    assert(ignoreFreq < 1. && 0. <= ignoreFreq); 
	  }
	  break; 
	case 'r':
	  {
	    constBurninWasSet  = true; 
	    auto &&iss = std::stringstream{optarg}; 
	    iss >> constBurnin ; 
	  }
	  break; 
	case 'm': 
	  takeAll = true; 
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
	    relBurninWasSet = true ; 
	    auto &&iss = std::istringstream {optarg}; 
	    iss >> relBurnin; 
	  }
	  break; 
	default : 
	  {
	    std::cerr << "unknown argument >" << c << "<" << std::endl; 
	    assert(0);
	  }
	} 
    }
  
  if(argc == 0 || files.size() == 0 )
    printUsage(); 

  for(auto file : files)
    {
      if(not std::ifstream(file)) 
	{
	  std::cout << "error: could not open file >" << file << "<" << std::endl; 
	  myExit(-1); 
	}
    }

  if(constBurninWasSet && relBurninWasSet)
    {
      std::cout << "Error: you cannot set the relative burn-in AND a range of trees to\n"
		<< "use.\n";
      myExit(-1); 
    }

  auto asdsf = SplitFreqAssessor(files); 
  
  nat end = asdsf.getMinNumTrees(); 
  
  if(end < constBurnin)
    {
      std::cout << "you are trying to discard " << constBurnin << " trees, but the minimum trees available in one of the files is just " << end << std::endl; 
      myExit(-1); 
    }

  nat start = 0;   
  if(constBurninWasSet)
    start = constBurnin; 
  else 
    start = nat(double(end ) * relBurnin); 

  asdsf.extractBipsNew(start , end , takeAll);
  auto asdsfResult = asdsf.computeAsdsfNew(ignoreFreq); 
  
  std::cout << "used trees number " << start << " to " <<  end << std::endl; 
  if(takeAll)
    std::cout << "but also included any tree with an id higher than " << end << ", if available." << std::endl; 
  std::cout << "average deviation of split frequencies: " <<  asdsfResult.first * 100  << "%" << std::endl; 
  std::cout  << "maximum deviation of split frequencies: " << asdsfResult.second * 100 << "%" << std::endl; 
  std::cout << "ignored splits that occurred in less than "  << ignoreFreq * 100 << "% of the trees for any of the specified files." << std::endl; 

// #ifdef _USE_GOOGLE_PROFILER
//   ProfilerStop();
// #endif

  return 0; 
}


