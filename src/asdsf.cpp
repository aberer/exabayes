#include "axml.h"
#include "proposalFunction.h"

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

// #undef PRINT 
// #define PRINT printf

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

  vector<string> tmp;
  for(int i = 3 ; i < argc; ++i)  
    {
      // check if the file exists 
      FILE *fh = fopen(argv[i],"r"); 
      if(fh == NULL)
	{
	  cerr << "could not find file >" << argv[i] << "<"<< endl; 
	  exit(0); 
	}

      tmp.push_back(string(argv[i]));
    }

  AvgSplitFreqAssessor asdsf(tmp); 

  asdsf.setStart(start);
  asdsf.setEnd(end);

  // cout << "all files contain at least " << asdsf.getMinNumTrees() << " trees" << endl; ;
  asdsf.extractBips();

  double ignoreFreq = 0.1; 
  
  cout << "average deviation of split frequencies: " << asdsf.computeAsdsf(ignoreFreq) * 100  << "%" << endl; 
  cout << "ignored splits that did not occur more than "  << ignoreFreq * 100 << "% of the trees for any of the specified files." << endl; 
  return 0; 
}

