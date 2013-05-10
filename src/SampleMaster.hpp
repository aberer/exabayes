/** 
    @file SampleMaster.hpp
    
    @brief represents various chains sampling the posterior probability space
    
    Despite of its modest name, this is in fact the master class.  
 */ 


#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <vector>

#include "axml.h"

#include "CommandLine.hpp"
#include "CoupledChains.hpp"
#include "ConfigReader.hpp"
#include "ParallelSetup.hpp"

using namespace std; 


class SampleMaster
{
public: 
  SampleMaster(const CommandLine &cl , ParallelSetup &pl);
  ~SampleMaster(){};

  void initRunParameters(string configFileName); 
  void finalizeRuns();  
  void run(); 
  void initWithConfigFile(string configFileName);
  void setupProposals(vector<AbstractProposal*> &result, const PriorBelief &priors ); 

  // HERE 
  void setGuidedRadius ( int guidedRadius ) {this->guidedRadius = guidedRadius ; }
  void setParsimonyWarp   ( double parsimonyWarp   ) {this->parsimonyWarp   = parsimonyWarp   ; }
  void setEsprStopProp ( int esprStopProp ) {this->esprStopProp = esprStopProp ; }
  void setTuneFreq  ( int tuneFreq  ) {this->tuneFreq  = tuneFreq  ; }

  void setNumGen(int numGen){this->numGen = numGen; }
  void setDiagFreq(int diagFreq) {this->diagFreq = diagFreq; }
  void setAsdsfIgoreFreq(double asdsfIgnoreFreq){this->asdsfIgnoreFreq = asdsfIgnoreFreq; }
  // void setMyBatch ( int myBatch ) {this->myBatch = myBatch ; }
  void setRunId ( string runId ) {this->runId = runId ; }
  void setNumRunConv ( int numRunConv ) {this->numRunConv = numRunConv ; }
  void setSamplingFreq ( int samplingFreq ) {this->samplingFreq = samplingFreq ; }
  void setBurninProportion ( double burninProportion ) {this->burninProportion = burninProportion ; }
  void setBurninGen ( int burninGen ) {this->burninGen = burninGen ; }
  void setAsdsfConvergence ( double asdsfConvergence ) {this->asdsfConvergence = asdsfConvergence ; }
  void setNumCoupledChains(int numCoupledChains) {this->numCoupledChains = numCoupledChains; }
  void setPrintFreq(int printFreq){this->printFreq = printFreq; }
  void setHeatFactor(double heat) {this->heatFactor = heat; }
  void setSwapInterval(int swapI){this->swapInterval = swapI; }
  void setTuneHeat(bool heat){this->tuneHeat = heat; }

  void validateRunParams(); 	// TODO  

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
  string runId; 
  int numCoupledChains; 
  int printFreq; 
  double heatFactor ; 
  int swapInterval; 
  bool tuneHeat; 
  int tuneFreq;  

  // move options
  int esprStopProp; 
  double parsimonyWarp;   
  int guidedRadius; 


  bool convergenceDiagnostic(); 
  // int myBatch; 

};  


#endif
