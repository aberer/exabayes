#include "axml.h"
#include "bayes.h"

#define _INCLUDE_DEFINITIONS
#include "globals.h"
#undef _INCLUDE_DEFINITIONS

#include "main-common.h"
#include "AvgSplitFreqAssessor.hpp"

#if HAVE_PLL != 0 
#include "globalVariables.h" 
#endif


int main(int argc, char** argv)
{
  if(argc < 4)
    {
      cout << "Usage: " << argv[0] << " start end [file ...]  " << endl << endl
	   << "where  " << endl
	   << "start\t is first tree to include (potentially skipping a burnin) " << endl
	   << "end\t is the last generation to include (file may be truncated)" << endl
	   << "[file ...]\t are various ExaBayes_topology* files. " << endl; 
      exit(0); 
    }

  int start = atoi(argv[1]); 
  int end = atoi(argv[2]); 
  
  char *aFile = argv[3]; 	// TODO 


  vector<string> tmp;
  tmp.push_back(string(aFile)); 
  

  AvgSplitFreqAssessor asdsf(tmp,start,end); 
  asdsf.extractBips();

  return 0; 
}
