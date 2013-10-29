#include <sstream>
#include <string.h>
#include <cassert>
#include <iostream>
#include <cmath>
#include "Branch.hpp"
#include "BasicTreeReader.hpp"
#include "TreeProcessor.hpp"
#include "parameters/BranchLengthsParameter.hpp"


TreeProcessor::TreeProcessor(std::vector<std::string> fileNames)  
{
  fillTaxaInfo(fileNames[0]); 

  initializeTreeOnly(taxa.size() );
  fns = fileNames; 
}


TreeProcessor::TreeProcessor(TreeProcessor&& tp) 
  : tralnPtr(std::move(tp.tralnPtr))
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


// static nodeptr getUnlinkedNode(const TreeAln &traln, nat id)
// {
//   auto p  = traln.getNode(id); 
//   if(p->back == NULL)
//     return p; 
//   else if(p->next->back == NULL)
//     return p->next; 
//   else if(p->next->next->back == NULL)
//     return p->next->next; 
//   else 
//     assert(0);
//   return NULL; 
// }



template<bool readBl>
void TreeProcessor::nextTree(std::istream &treefile) 
{
  auto paramPtr = std::unique_ptr<AbstractParameter>(new BranchLengthsParameter(0,0));   
  paramPtr->addPartition(0);

  while( treefile.get() != '('); 
  treefile.unget();

  auto bt = BasicTreeReader<IntegerLabelReader,
			    typename std::conditional<readBl,
						      ReadBranchLength,
						      IgnoreBranchLength>::type>(taxa.size());

  auto branches = bt.extractBranches(treefile);

  tralnPtr->unlinkTree();
  for(auto b :branches)
    {
      tralnPtr->clipNode(tralnPtr->getUnhookedNode(b.getPrimNode()), tralnPtr->getUnhookedNode(b.getSecNode()) );
      if(readBl)
      	tralnPtr->setBranch(b, paramPtr.get());
    }
}

void TreeProcessor::skipTree(std::istream &iss)
{
  while( iss &&  iss.get() != ';'); 
}


void TreeProcessor::initializeTreeOnly(int numTax )
{
  tralnPtr = std::unique_ptr<TreeAln>(new TreeAln());
  tree *tr = tralnPtr->getTr();
  tr->mxtips = numTax; 

#if HAVE_PLL != 0
  partitionList *pl = (partitionList*)exa_calloc(1,sizeof(partitionList)); 
  // pl->numberOfPartitions = 1; 	// BAD!!!! needed for extractBips
  setupTree(tr, false, pl);  
  tralnPtr->setPartitionList(pl); 
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



template void TreeProcessor::nextTree<true>(std::istream &treefile) ; 
template void TreeProcessor::nextTree<false>(std::istream &treefile) ;
