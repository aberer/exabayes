/**
   @file PriorManager.hpp

   @brief top level object that represents the prior probability of a
   chain
 */ 


#ifndef _PRIORMANAGER_H
#define _PRIORMANAGER_H

#include "nclConfigReader.h"

using namespace std; 

#include <vector>





// okay, since we have so many priors and they (mostly) do so few
// things, let's all keep them here: 

class BranchPrior
{
  
}; 


class RevMatPrior 
{
}; 


class FreqPrior
{
  
}; 


class RateHetPrior
{
  
}; 


// ################################################################


class PriorManager
{
public: 
  PriorManager(initParamStruct &initParams); 
  ~PriorManager(){};

  // TODO support topology 

  void updateBranchLength(double oldValue, double newValue); 
  void updateRevMat( vector<double> oldValues, vector<double> newValues); 
  void updateFreq(vector<double> oldValues, vector<double> newValues); 
  void updateRateHet(double oldValue, double newValue); 
  
  double getLogProb(){return logProb;}
  
private: 
  double logProb; 
  
  BranchPrior bP;
  RevMatPrior rvP; 
  FreqPrior fP; 
  RateHetPrior rhP; 

}; 


#endif
