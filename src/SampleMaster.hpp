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

using namespace std; 


class SampleMaster
{
public: 
  SampleMaster(const ParallelSetup &pl, const CommandLine& cl ) ; 

  void initializeRuns( ); 
  void initRunParameters(string configFileName); 
  void finalizeRuns();  
  void run(); 
  void initWithConfigFile(string configFileName, const TreeAlnPtr &traln, vector<unique_ptr<AbstractProposal> > &proposalResult, vector<shared_ptr<RandomVariable> > &variableResult, LikelihoodEvaluatorPtr &eval); 
  void validateRunParams(); 	// TODO  
  void branchLengthsIntegration()  ;  

private: 
  void initTrees(vector<TreeAlnPtr> &trees ); // METHODS 


private:			// ATTRIBUTES 
  bool convergenceDiagnostic(); 

  vector<CoupledChainsPtr> runs; // TODO bad design: just want to avoid getting memory leaks

  ParallelSetup pl; 

  CLOCK::system_clock::time_point initTime; 

  BlockParams paramBlock; 
  BlockRunParameters runParams;  
  BlockProposalConfig propConfig;   

  Randomness masterRand; 
  const CommandLine cl; 
};  

#endif
