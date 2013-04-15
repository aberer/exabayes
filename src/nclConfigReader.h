
/**
   @file nclConfigReader.h
   @brief Wraps the ncl-nexus parser for usage with ExaBayes.     
*/

#ifndef _NCL_CONFIG_READER
#define _NCL_CONFIG_READER

#include "proposalType.h"

typedef struct  
{
  double initWeights[NUM_PROPOSALS]; 
  double initPenaltyFactor;   
  double eSprStopProb; 
  int numGen; 
  int samplingFrequency; 
  int numIndiChains; 
  int diagFreq; 
  int numCoupledChains; 
  int initGuidedSPR; 
  int printFreq; 
  double asdsfIgnoreFreq; 
  double asdsfConvergence; 
  double heatFactor; 
  int swapInterval; 
  bool tuneHeat;
  int burninGen ; 
  double burninProportion; 
  int tuneFreq; 
} initParamStruct ; 


#ifdef __cplusplus
extern "C"{
#endif
  void parseConfigWithNcl(char *configFileName, initParamStruct **initParam); 
#ifdef __cplusplus
}
#endif

#endif
