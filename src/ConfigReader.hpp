
/**
   @file nclConfigReader.h
   @brief Wraps the ncl-nexus parser for usage with ExaBayes.     
*/

#ifndef _NCL_CONFIG_READER
#define _NCL_CONFIG_READER

#include <ncl/ncl.h>
#include <vector>

#include "proposalType.h"
#include "PriorManager.hpp"
#include "AbstractProposal.hpp"

class ConfigReader : public NxsReader
{
public: 
  ConfigReader() : NxsReader(){SetWarningOutputLevel(SUPPRESS_WARNINGS_LEVEL); }
  virtual void ExitingBlock(NxsString blockName){}
  virtual void ExecuteStopping(){}
  virtual void ExecuteStarting(){}
  
  PriorManager getPriorConfig(){return priorConfig; }
  

private: 
  PriorManager priorConfig;
  vector<AbstractProposal*>  proposals;     
  int numCoupledChains; 
}; 


typedef struct  
{
  double initWeights[NUM_PROPOSALS]; 
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
  int numRunParallel; 
  double parsWarp; 
} initParamStruct ; 

void parseConfigWithNcl(char *configFileName, initParamStruct **initParam); 

#endif
