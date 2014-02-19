/** 
    @file SplitFreqAssessor.hpp

    @brief Calculates the asdsf. 

    @notice This file is full of hacks, to get this somehow going. 
 */ 

#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

// #include "axml.h"

#include "tree-parse/TreeProcessor.hpp"
#include "BipartitionHash.hpp"

class SplitFreqAssessor : public TreeProcessor 
{
public: 
  SplitFreqAssessor(std::vector<std::string>fileNames, bool expensiveCheck);
  /** 
      @brief use the new bipartition hash for extracting bipartitions 
   */ 
  void extractBipsNew(nat start, nat end, bool takeAll); 
  /** 
      @brief gets the minimum number of trees present in all of the files 
  */ 
  int getMinNumTrees(); 
  
  std::pair<double,double> computeAsdsfNew(double ignoreFreq);

private: 			// ATTRIBUTES
  nat getNumTreeAvailable(std::string filename); 
  std::vector<BipartitionHash> newBipHashes;   
  std::unordered_map<std::string, nat> file2numTree; 
}; 

#endif
