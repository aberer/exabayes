#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

#include "axml.h"
#include "BipartitionHash.hpp"

/**
   @brief a one-time object that computes the ASDSF for the trees
   contained in a bunch of files.

 */ 
class AvgSplitFreqAssessor
{
public: 
  AvgSplitFreqAssessor(vector<string>fileNames);
  ~AvgSplitFreqAssessor();

  /** @brief return the asdsf of the respective trees in the respective range */ 
  double computeAsdsf(double ignoreFreq); 
  
  /** @brief add bipartitions in the current traln structure into the bipartition hash */
  void extractBips(); 

  /** @brief gets the minimum number of trees present in all of the files */ 
  int getMinNumTrees(); 

  int getEnd(){return end; } 
  int getStart(){return start; }
  void setEnd(int _end){end = _end; }
  void setStart(int _start){start = _start; }

  int getNumTreeAvailable(string filename); 

  static double relativeWeight;

private: 
  void fillTaxaInfo(string fileName); 
  bool fileIsCorrect(string fileName);   
  void nextTree(FILE *fh);
  void initializeTreeOnly(int numTax); 

  TreeAln *traln;
  vector<string> fns; 
  vector<string>  taxa; 
  int start; 
  int end; 
  
  BipartitionHash* bipHash;
}; 

#endif
