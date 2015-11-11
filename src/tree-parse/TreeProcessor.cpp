#include <sstream>
#include "extensions.hpp" 
#include <string.h>
#include <cassert>
#include "TreeInitializer.hpp"
#include <iostream>
#include <cmath>
#include "Branch.hpp"
#include "BasicTreeReader.hpp"
#include "TreeProcessor.hpp"
#include "parameters/BranchLengthsParameter.hpp"


TreeProcessor::TreeProcessor(std::vector<std::string> fileNames)  
{
  assert(fileNames.size() > 0); 
  fillTaxaInfo(fileNames.at(0)); 
  nat numTax = taxa.size(); 
  tralnPtr = std::unique_ptr<TreeAln>(new TreeAln(numTax));
  TreeInitializer::initializeBranchLengths(tralnPtr->getTrHandle(), 1,numTax); 
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


template<bool readBl>
void TreeProcessor::nextTree(std::istream &treefile) 
{
  auto paramPtr = make_unique<BranchLengthsParameter>(0,0, std::vector<nat>{0});   
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
	{
	  tralnPtr->setBranchUnchecked(b);
	}
    }
}

void TreeProcessor::skipTree(std::istream &iss)
{
  while( iss &&  iss.get() != ';'); 
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
