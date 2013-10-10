/** 
    @file SplitFreqAssessor.hpp

    @brief Calculates the asdsf. 

    @notice This file is full of hacks, to get this somehow going. 
 */ 

#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

#include "axml.h"

#include "TreeProcessor.hpp"
class BipartitionHashNew; 

class SplitFreqAssessor : public TreeProcessor 
{
public: 
  SplitFreqAssessor(std::vector<std::string>fileNames);
  ~SplitFreqAssessor();
  /** 
      @brief use the new bipartition hash for extracting bipartitions 
   */ 
  void extractBipsNew(); 
  /** 
      @brief gets the minimum number of trees present in all of the files 
  */ 
  int getMinNumTrees(); 
  
  std::pair<double,double> computeAsdsfNew(double ignoreFreq);
  
  int getEnd(){return end; } 
  int getStart(){return start; }
  void setEnd(int _end){end = _end; }
  void setStart(int _start){start = _start; }

  int getNumTreeAvailable(std::string filename); 

private: 			// ATTRIBUTES
  int start; 
  int end; 
  std::vector<BipartitionHashNew> newBipHashes;   
}; 

#endif
