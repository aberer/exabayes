/** 
    @file SampleMaster.hpp
    
    @brief represents various chains sampling the posterior probability space
    
    Despite of its modest name, this is in fact the master class.  
 */ 

#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <vector>

#include "axml.h"
#include "config/BlockRunParameters.hpp"
#include "config/BlockProposalConfig.hpp"
#include "config/CommandLine.hpp"
#include "CoupledChains.hpp"
#include "config/ConfigReader.hpp"
#include "ParallelSetup.hpp"
#include "time.hpp"
#include "Checkpointable.hpp"
#include "DiagnosticsFile.hpp"


class SampleMaster : public Checkpointable
{
public: 
  SampleMaster(const ParallelSetup &pl, const CommandLine& cl ) ; 

  void initializeRuns( ); 
  void initRunParameters(string configFileName); 
  void finalizeRuns();  
  void run(); 
  void initWithConfigFile(string configFileName, shared_ptr<TreeAln> traln, 
			  vector<unique_ptr<AbstractProposal> > &proposalResult, vector<unique_ptr<AbstractParameter> > &variableResult, 
			  shared_ptr<LikelihoodEvaluator> eval); 
  void validateRunParams(); 	// TODO  
  void branchLengthsIntegration()  ;  

  void printAlignmentInfo(const TreeAln &traln); 

  virtual void readFromCheckpoint( std::ifstream &in ) ; 
  virtual void writeToCheckpoint( std::ofstream &out)  ;   

private: 
  void writeCheckpointMaster(); 
  double convergenceDiagnostic(); 
  void initTrees(vector<shared_ptr<TreeAln> > &trees, randCtr_t seed, nat &treesConsumed, nat numTreesAvailable, FILE *fh); 
  void printDuringRun(nat gen); 

private:			// ATTRIBUTES 
  vector<CoupledChains> runs; 
  ParallelSetup pl; 
  CLOCK::system_clock::time_point initTime; 
  BlockParams paramBlock; 
  BlockRunParameters runParams;  
  BlockProposalConfig propConfig;   
  Randomness masterRand;   	// not checkpointed
  CommandLine cl; 
  CLOCK::system_clock::time_point timeIncrement; 
  DiagnosticsFile diagFile; 
};  

#endif
