
/** 
    @file AvgSplitFreqAssessor.hpp

    @brief Calculates the asdsf. 

    @notice This file is full of hacks, to get this somehow going. 
 */ 

#include <fstream>
#include <algorithm>
#include <sstream> 

#include <ncl/ncl.h>

#include "axml.h"
#include "bayes.h"
#include "TreeAln.hpp"

#include "output.h"
#include "treeRead.h"

#include "AvgSplitFreqAssessor.hpp"



/**
   @brief checks, if everything is in order with the files and the trees it contains  
*/ 
AvgSplitFreqAssessor::AvgSplitFreqAssessor(vector<string> fileNames)
{
  fillTaxaInfo(fileNames[0]); 

  // at least does some checking 
  for(string fn : fileNames)
    getNumTreeAvailable (fn); 

  initializeTreeOnly(taxa.size() );
  fns = fileNames; 

  bipHash = new BipartitionHash(taxa.size(), fileNames.size());

  this->start = 0; 
  this->end = getMinNumTrees();
  // cout << "AvgSplitFreqAssessor: start " << start << "\tend" <<  end << endl; 
}




AvgSplitFreqAssessor::~AvgSplitFreqAssessor()
{
  delete bipHash; 
  delete traln; 
}


void AvgSplitFreqAssessor::nextTree(FILE *fh )
{
  tree *tr = traln->getTr();
  int c = 0; 
  while( (c = getc(fh)) != '('); 
  ungetc(c, fh);   
  myTreeReadLen(fh, tr , FALSE); 
}




void AvgSplitFreqAssessor::extractBips()
{
  int ctr = 0; 
  for (auto filename : fns)
    {
      FILE *fh = fopen(filename.c_str(), "r"); 

      for(int i = 0; i < start; ++i)
	nextTree(fh);

      for(int i = start ; i < end; ++i)
	{
	  nextTree(fh);
	  bipHash->addBipartitionsToHash(*traln, ctr);      
	}
      
      fclose(fh);
      ++ctr; 
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
// bool AvgSplitFreqAssessor::fileIsCorrect(string fileName)
// {
//   int numTrees = getNumTreeAvailable(fileName);
//   return (numTrees >=  (end - start)); 
// }


double AvgSplitFreqAssessor::computeAsdsf(double ignoreFreq)
{
  return bipHash->averageDeviationOfSplitFrequencies(ignoreFreq); 
}




int AvgSplitFreqAssessor::getNumTreeAvailable(string fileName)
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
  return numTrees; 
}



void AvgSplitFreqAssessor::fillTaxaInfo(string fileName)
{
  string whiteSpace = " \t"; 
  // cout << "filename is " << fileName << endl;
  ifstream infile(fileName); 
  string line; 
  bool foundStart = false; 
  bool abort = false; 
  while(not abort && getline(infile, line))
    {      
      string cleanLine = trim(line); 

      // cout << "line is  >" << line << "<" << endl; 

      if(foundStart)
	{
	  if(cleanLine[cleanLine.size()-1] == ';') // we are done 
	    abort = true; 

	  string cleanerString = trim(cleanLine, ",;"); 

	  int pos = cleanerString.find_first_of(whiteSpace, 0);
	  string num =  cleanerString.substr(0, pos),
	    name = cleanerString.substr(pos+1, cleanerString.size()); 	  

	  taxa.push_back(name); 
	}
      else if(cleanLine.compare("translate") == 0  )
	foundStart = true; 
    }  

  assert(foundStart); 
}



/** 
    @brief only initializes a raw tree, no partitions or alignment information. 
    
    important: does NOT need a bytefile 
 */ 
void AvgSplitFreqAssessor::initializeTreeOnly(int numTax )
{
  traln = new TreeAln();
  tree *tr = traln->getTr();
  tr->mxtips = numTax; 

#if HAVE_PLL != 0
  partitionList *pl = (partitionList*)exa_calloc(1,sizeof(partitionList)); 
  pl->numberOfPartitions = 0; 
  setupTree(tr, false, pl);  
  traln->setPartitionList(pl); 
#else 
  tr->NumberOfModels = 0; 
  setupTree(tr);
#endif

  int space = int(log(numTax) * 10 ) ;
  // initialize the name hash, s.t. we can read trees 
  for(int i = 1; i <= tr->mxtips; i++)
    {      
      tr->nameList[i] = (char*)malloc(sizeof(char) * space );      
      stringstream ss; 
      ss << i ;
      strcpy(tr->nameList[i], ss.str().c_str()); 
      addword(tr->nameList[i] , tr->nameHash, i); 
    }
}




int AvgSplitFreqAssessor::getMinNumTrees()
{
  int minimum = -1; 
  for(auto fn : fns )
    {
      FILE *fh = fopen(fn.c_str(), "r"); 

      // the following is moderately evil: assuming, there is no '('
      // except for the first tree and no ';' char except for tree
      // delimitation, after we found a 

      int numTreesHere = 0; 

      int c; 
      while( ( c = getc(fh) ) != EOF && c != '('); 
      
      
      while( (c = getc(fh) )  != EOF )
	if(c == ';')
	  numTreesHere++;

      if(minimum == -1 || numTreesHere < minimum)
	minimum = numTreesHere; 

      fclose(fh); 
    }
  
  return minimum; 
} 
