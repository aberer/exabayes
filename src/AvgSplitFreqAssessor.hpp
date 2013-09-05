/** 
    @file AvgSplitFreqAssessor.hpp

    @brief Calculates the asdsf. 

    @notice This file is full of hacks, to get this somehow going. 
 */ 

#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

#include "axml.h"
// #include "BipartitionHash.hpp"

#include "TreeProcessor.hpp"
class BipartitionHashNew; 

class AvgSplitFreqAssessor : public TreeProcessor 
{
public: 
  AvgSplitFreqAssessor(std::vector<std::string>fileNames);
  ~AvgSplitFreqAssessor();
  /** 
      @brief return the asdsf of the respective trees in the respective range 
  */ 
  // double computeAsdsf(double ignoreFreq); 
  /** 
      @brief add bipartitions in the current traln structure into the
      bipartition hash
  */
  // void extractBips(); 
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
