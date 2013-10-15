#include <sstream>
#include <string.h>
#include <cassert>
#include <iostream>
#include <cmath>

#include "TreeProcessor.hpp"
#include "treeRead.h"

TreeProcessor::TreeProcessor(std::vector<std::string> fileNames)  
{
  fillTaxaInfo(fileNames[0]); 

  initializeTreeOnly(taxa.size() );
  fns = fileNames; 
}


TreeProcessor::TreeProcessor(TreeProcessor&& tp) 
  : traln(std::move(tp.traln))
  , fns(tp.fns)
  , taxa(tp.taxa)
{
  
}  


TreeProcessor& TreeProcessor::operator=(TreeProcessor &&rhs)
{
  if(this == &rhs)
    return *this; 
  else 
    assert(0); 
} 


void TreeProcessor::nextTree(std::istream &treefile, bool readBL) 
{
  tree *tr = traln->getTr();
  while( treefile.get() != '('); 
  treefile.unget();
  treefile.unget();

  auto treestring = std::string {}; 
  std::getline(treefile, treestring); 
  // std::cout << "reading " << treestring << std::endl; 
  myTreeReadLen(treestring, tr , readBL ? TRUE : FALSE); 
}


void TreeProcessor::initializeTreeOnly(int numTax )
{
  traln = std::unique_ptr<TreeAln>(new TreeAln());
  tree *tr = traln->getTr();
  tr->mxtips = numTax; 

#if HAVE_PLL != 0
  partitionList *pl = (partitionList*)exa_calloc(1,sizeof(partitionList)); 
  // pl->numberOfPartitions = 1; 	// BAD!!!! needed for extractBips
  setupTree(tr, false, pl);  
  traln->setPartitionList(pl); 
#else 
  tr->NumberOfModels = 0; 
  setupTree(tr);
#endif

  for(int i = 0; i < numTax + 3 * ( numTax - 1  ) ; ++i ) 
    {
      tr->nodeBaseAddress[i].z = (double*)exa_malloc(1 * sizeof(double)); 
    }

  int space = int(log(numTax) * 10 ) ;
  // initialize the name hash, s.t. we can read trees 
  for(int i = 1; i <= tr->mxtips; i++)
    {      
      tr->nameList[i] = (char*)malloc(sizeof(char) * space );      
      auto &&ss = std::stringstream{}; 
      ss << i ;
      strcpy(tr->nameList[i], ss.str().c_str()); 
      addword(tr->nameList[i] , tr->nameHash, i); 
    }
}


std::string TreeProcessor::trim(const std::string& str, const std::string& whitespace) 
{
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; 
  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}


void TreeProcessor::fillTaxaInfo(std::string fileName)
{
  std::string whiteSpace = " \t"; 
  std::ifstream infile(fileName); 
  std::string line; 
  bool foundStart = false; 
  bool abort = false; 
  while(not abort && getline(infile, line))
    {      
      std::string cleanLine = trim(line); 

      if(foundStart)
	{
	  if(cleanLine[cleanLine.size()-1] == ';') // we are done 
	    abort = true; 

	  std::string cleanerString = trim(cleanLine, ",;"); 

	  int pos = cleanerString.find_first_of(whiteSpace, 0);
	  std::string num =  cleanerString.substr(0, pos),
	    name = cleanerString.substr(pos+1, cleanerString.size()); 	  

	  taxa.push_back(name); 
	}
      else if(cleanLine.compare("translate") == 0  )
	foundStart = true; 
    }  

  assert(foundStart); 
}
