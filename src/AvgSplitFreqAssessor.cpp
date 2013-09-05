
#include <fstream>
#include <algorithm>
#include <sstream> 
#include <cassert>

#include <ncl/ncl.h>

#include "axml.h"
#include "TreeAln.hpp"
#include "Arithmetics.hpp"
#include "AvgSplitFreqAssessor.hpp"
#include "BipartitionHash.hpp"


AvgSplitFreqAssessor::AvgSplitFreqAssessor(vector<string> fileNames)
  : TreeProcessor(fileNames)
{
  for(string fn : fileNames)
    getNumTreeAvailable (fn); 

  this->start = 0; 
  this->end = getMinNumTrees();

  // bipHash = new BipartitionHash(taxa.size(), fileNames.size());
}


AvgSplitFreqAssessor::~AvgSplitFreqAssessor()
{
  // delete bipHash; 
}


void AvgSplitFreqAssessor::extractBipsNew()
{
  int ctr = 0; 
  for (auto filename : fns)
    {
      FILE *fh = fopen(filename.c_str(), "r"); 

      for(int i = 0; i < start; ++i)
	nextTree(fh);

      auto bipHash = BipartitionHashNew(traln->getNumberOfTaxa()); 
      for(int i = start ; i < end; ++i)
	{
	  nextTree(fh);
	  bipHash.addTree(*traln,false);
	}
      newBipHashes.push_back(bipHash);
      
      fclose(fh);
      ++ctr; 
    }
}

// void AvgSplitFreqAssessor::extractBips()
// {
//   int ctr = 0; 
//   for (auto filename : fns)
//     {
//       FILE *fh = fopen(filename.c_str(), "r"); 

//       for(int i = 0; i < start; ++i)
// 	nextTree(fh);

//       for(int i = start ; i < end; ++i)
// 	{
// 	  nextTree(fh);
// 	  // bipHash->addBipartitionsToHash(*traln, ctr);      
// 	}
      
//       fclose(fh);
//       ++ctr; 
//     }
// }


// double AvgSplitFreqAssessor::computeAsdsf(double ignoreFreq)
// {
//   return bipHash->averageDeviationOfSplitFrequencies(ignoreFreq); 
// }

auto AvgSplitFreqAssessor::computeAsdsfNew(double ignoreFreq)
 -> std::pair<double,double> 
{
  auto allBips = std::unordered_set<Bipartition>(); 

  nat numTrees = newBipHashes[0].getTreesAdded();
  for(auto &elem : newBipHashes)
    assert(numTrees == elem.getTreesAdded()); 
  
  for(auto &h : newBipHashes)
    for(auto &elem : h)
      allBips.insert(elem.first);

  std::unordered_map<Bipartition,std::vector<double> > numOccs; 

  for(auto &elem : allBips)
    {
      auto tmp = std::vector<double>{};
      for(auto &h : newBipHashes)
	tmp.push_back( double(h.getPresence(elem).count()) );
      numOccs[elem] = tmp; 
    }
  
  auto asdsfPart = std::vector<double>{};
  for(auto &elem : allBips)
    {
      bool isRelevant = false;  
      for(auto &v : numOccs[elem])
	{
	  v /= numTrees; 
	  isRelevant |= (ignoreFreq < v); 
	}

      auto var = Arithmetics::getVariance(numOccs[elem]); 

      if(isRelevant)
	asdsfPart.push_back(sqrt(var));
    }

  auto asdsfMax = std::max_element(asdsfPart.begin(), asdsfPart.end()); 
  auto mean = Arithmetics::getMean(asdsfPart); 

  return make_pair(mean, *asdsfMax); 
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
      std::string cleanline = this->trim(line); 

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
