#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

#include "axml.h"
#include "BipartitionHash.hpp"
#include "BipartitionHashNew.hpp"

/**
   @brief a one-time object that computes the ASDSF for the trees
   contained in a bunch of files.

 */ 
class AvgSplitFreqAssessor
{
public: 
  AvgSplitFreqAssessor(std::vector<std::string>fileNames);
  ~AvgSplitFreqAssessor();
  /** 
      @brief return the asdsf of the respective trees in the respective range 
  */ 
  double computeAsdsf(double ignoreFreq); 
  /** 
      @brief add bipartitions in the current traln structure into the
      bipartition hash
  */
  void extractBips(); 
  /** 
      @brief use the new bipartition hash for extracting bipartitions 
   */ 
  void extractBipsNew(); 
  /** 
      @brief gets the minimum number of trees present in all of the files 
  */ 
  int getMinNumTrees(); 
  
  double computeAsdsfNew(double ignoreFreq);

  int getEnd(){return end; } 
  int getStart(){return start; }
  void setEnd(int _end){end = _end; }
  void setStart(int _start){start = _start; }

  int getNumTreeAvailable(std::string filename); 

  static double relativeWeight;

private: 
  void fillTaxaInfo(std::string fileName); 
  bool fileIsCorrect(std::string fileName);   
  void nextTree(FILE *fh);
  void initializeTreeOnly(int numTax); 

  TreeAln *traln;
  std::vector<std::string> fns; 
  std::vector<std::string>  taxa; 
  int start; 
  int end; 
  
  BipartitionHash* bipHash;

  std::vector<BipartitionHashNew> newBipHashes;   
}; 

#endif
