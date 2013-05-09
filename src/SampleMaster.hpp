/** 
    @file SampleMaster.hpp
    
    @brief represents various chains sampling the posterior probability space
    
    Despite of its modest name, this is in fact the master class.  
 */ 


#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <vector>

// #include "bayes.h"
#include "axml.h"

#include "CoupledChains.hpp"

using namespace std; 


class SampleMaster
{
public: 
  SampleMaster(analdef *adef , int seed, initParamStruct *initParams );
  ~SampleMaster(){};

  void finalizeRuns();  
  void  run(); 

private: 
  vector<CoupledChains> runs; 
  int diagFreq ; 
  double asdsfIgnoreFreq; 
  double asdsfConvergence; 
  int burninGen; 
  double burninProportion; 
  int samplingFreq; 
  int numRunConv; 
  int numGen; 

  bool convergenceDiagnostic(); 

};  


#endif
