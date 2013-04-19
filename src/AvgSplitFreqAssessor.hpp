#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

#include "axml.h"
#include "bayes.h"


#include <vector>



/**
   @brief a one-time object that computes the ASDSF for the trees
   contained in a bunch of files.

 */ 
class AvgSplitFreqAssessor
{
public: 
  AvgSplitFreqAssessor(vector<string>fileNames,  int start, int end);
  double computeAsdsf(); 
  void extractBips(); 

private: 
  void fillTaxaInfo(string fileName); 
  bool fileIsCorrect(string fileName);   

  tree *tr; 
  vector<string> fns; 
  vector<string>  taxa; 
  int start; 
  int end; 
  


}; 

#endif
