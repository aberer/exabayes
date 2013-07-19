#include <cassert>

#include "axml.h"

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "AvgSplitFreqAssessor.hpp"

#if HAVE_PLL != 0 
#include "globalVariables.h" 
#endif


void genericExit(int code); 


int main(int argc, char** argv)
{
  if(argc < 4)
    {
      std::cout << "Usage: " << argv[0] << " start end [file ...]  " << std::endl << std::endl
	   << "where  " << std::endl
	   << "start\t is first tree to include (potentially skipping a burnin) " << std::endl
	   << "end\t is the last generation to include (file may be truncated)" << std::endl
	   << "[file ...]\t are various ExaBayes_topology* files. " << std::endl; 
      exit(-1); 
    }

  int start = atoi(argv[1]); 
  int end = atoi(argv[2]); 

  std::vector<std::string> tmp;
  for(int i = 3 ; i < argc; ++i)  
    {
      // check if the file exists 
      FILE *fh = fopen(argv[i],"r"); 
      if(fh == NULL)
	{
	  std::cerr << "could not find file >" << argv[i] << "<"<< std::endl; 
	  exit(-1); 
	}
      
      tmp.push_back(std::string(argv[i]));
    }

  AvgSplitFreqAssessor asdsf(tmp); 

  asdsf.setStart(start);
  asdsf.setEnd(end);

  // cout << "all files contain at least " << asdsf.getMinNumTrees() << " trees" << endl; ;
  asdsf.extractBips();

  double ignoreFreq = 0.1; 
  
  std::cout << "average deviation of split frequencies: " << asdsf.computeAsdsf(ignoreFreq) * 100  << "%" << std::endl; 
  std::cout << "ignored splits that did not occur more than "  << ignoreFreq * 100 << "% of the trees for any of the specified files." << std::endl; 
  return 0; 
}


