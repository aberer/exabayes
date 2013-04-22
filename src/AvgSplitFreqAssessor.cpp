
/** 
    @file AvgSplitFreqAssessor.hpp

    @brief Calculates the asdsf. 

    @notice This file is full of hacks, to get this somehow going. 
 */ 


#include <fstream>
#include <ncl/ncl.h>
#include <algorithm> 

#include "axml.h"
#include "bayes.h"
#include "main-common.h"

#include "AvgSplitFreqAssessor.hpp"


static void initializeTreeOnly(int numTax, tree **tre ); 


/**
   @brief checks, if everything is in order with the files and the trees it contains  
*/ 
AvgSplitFreqAssessor::AvgSplitFreqAssessor(vector<string> fileNames, int _start, int _end)
  : start(_start)
  , end(_end)
{
  fillTaxaInfo(fileNames[0]); 

  for(string fn : fileNames)
    {
      if(not fileIsCorrect(fn))
	{
	  cerr << "problem parsing file " << fn << endl;  
	  assert(0); 
	}
    }

  tr = NULL; 
  initializeTreeOnly(taxa.size(), &tr );
  unsigned int vectorLength = 0; 
  tr->bitVectors = initBitVector(tr->mxtips, &vectorLength );
  tr->h = initHashTable(10 * tr->mxtips);
  fns = fileNames; 
}




#include "output.h"
#include "treeRead.h"


void AvgSplitFreqAssessor::extractBips()
{
  for (auto filename : fns)
    {
      FILE *fh = fopen(filename.c_str(), "r"); 
      int c = 0; 
      while( (c = getc(fh)) != '('); 
      ungetc(c, fh); 

      myTreeReadLen(fh, this->tr, FALSE); 
      Tree2stringNexus(tr->tree_string, tr,  tr->nodep[1]->back, 0 ); 
      cout << tr->tree_string << endl; 

      exit(0); 

    }
}





std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; 
  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}



/**
   @brief checks, if the number <-> taxon mapping is correct and
   whether there are as many trees as requested for the asdsf
   
   @notice I guess this is not how you write a good parser ;-) 
 */
bool AvgSplitFreqAssessor::fileIsCorrect(string fileName)
{
  string whitespace = " \t"; 
  ifstream infile(fileName); 
  string line ; 
  bool foundTaxaStart = false,
    foundTreeStart = false; 
  int numTrees = 0; 

  while(getline(infile, line))
    {
      string cleanline = trim(line); 

      if(foundTreeStart)	// check number of trees 
	{
	  if(cleanline[cleanline.size()-1] == ';') 
	    ++numTrees; 
	  else 
	    break; 		// DONE  
	}
      else if(foundTaxaStart) 	// check if all taxa are there and have the appropriate number 
	{

	  if(cleanline[cleanline.size()-1] == ';')
	    foundTreeStart = true ; 

	  string cleanerString =  trim(cleanline,  ",;"); 
	  int pos = cleanerString.find_first_of(whitespace, 0); 
	  string num = cleanerString.substr(0, pos),
	    name = cleanerString.substr(pos+1, cleanerString.size()); 

	  int index; 
	  istringstream(num) >> index ; 
	  index--; 
	  assert(name.compare(taxa[index]) == 0); // assert same taxa names 
	  
	} 
      else if(cleanline.compare("translate") == 0)
	foundTaxaStart = true; 
    }

  assert(foundTaxaStart && foundTreeStart); 
  return (numTrees >=  (end - start)); 
}



double AvgSplitFreqAssessor::computeAsdsf()
{
  return 0; 
}


void AvgSplitFreqAssessor::fillTaxaInfo(string fileName)
{
  string whiteSpace = " \t"; 
  ifstream infile(fileName); 
  string line; 
  bool foundStart = false; 
  while(getline(infile, line))
    {      
      string cleanLine = trim(line); 

      if(foundStart)
	{
	  bool abort = false; 
	  if(cleanLine[cleanLine.size()-1] == ';') // we are done 
	    abort = true; 

	  string cleanerString = trim(cleanLine, ",;"); 

	  int pos = cleanerString.find_first_of(whiteSpace, 0);
	  string num =  cleanerString.substr(0, pos),
	    name = cleanerString.substr(pos+1, cleanerString.size()); 	  

	  taxa.push_back(name); 
	  if(abort)
	    break; 
	}
      else if(  cleanLine.compare("translate") == 0  )
	foundStart = true; 
    }  
  assert(foundStart); 
}


/** 
    @brief only initializes a raw tree, no partitions or alignment information. 
    
    important: does NOT need a bytefile 
 */ 
static void initializeTreeOnly(int numTax, tree **tre )
{
  *tre = (tree*)exa_calloc(1,sizeof(tree)); 
  tree *tr = *tre; 
  
  tr->mxtips = numTax; 

#if HAVE_PLL != 0
  partitionList pl; 
  
  setupTree(tr, false, &pl);
#else 
  setupTree(tr);
#endif

}

