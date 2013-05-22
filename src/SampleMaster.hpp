/** 
    @file SampleMaster.hpp
    
    @brief represents various chains sampling the posterior probability space
    
    Despite of its modest name, this is in fact the master class.  
 */ 


#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <vector>

#include "axml.h"
#include "BlockRunParameters.hpp"
#include "BlockProposalConfig.hpp"
#include "CommandLine.hpp"
#include "CoupledChains.hpp"
#include "ConfigReader.hpp"
#include "ParallelSetup.hpp"

using namespace std; 


class SampleMaster
{
public: 
  SampleMaster( const ParallelSetup &pl) ;   
  ~SampleMaster(){};

  void initializeRuns(const CommandLine &cl ); 
  void initRunParameters(string configFileName); 
  void finalizeRuns();  
  void run(); 
  void initWithConfigFile(string configFileName, const TreeAln &traln, vector<unique_ptr<AbstractProposal> > &proposalResult, vector<RandomVariable> &variableResult); 
  void validateRunParams(); 	// TODO  

private: 
  void initTrees(vector<TreeAln*> &trees, const CommandLine &cl ); 
  bool convergenceDiagnostic(); 

  vector<CoupledChains> runs; 

  ParallelSetup pl; 

  double initTime; 

  BlockRunParameters runParams;  
  BlockProposalConfig propConfig; 
  
};  

#endif
