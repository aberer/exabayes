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
#include "BlockRunParameters.hpp"

using namespace std; 


class SampleMaster
{
public: 
  SampleMaster(const CommandLine &cl , const ParallelSetup &pl) ;   
  ~SampleMaster(){};

  void initRunParameters(string configFileName); 
  void finalizeRuns();  
  void run(); 
  void initWithConfigFile(string configFileName, PriorBelief &prior, vector<double> &proposalWeights, const TreeAln& traln ); 
  void setupProposals(vector<Category> &proposalCategories, vector<double> proposalWeights, const PriorBelief &prior);

  // HERE 
  void setGuidedRadius ( int guidedRadius ) {this->guidedRadius = guidedRadius ; }
  void setParsimonyWarp   ( double parsimonyWarp   ) {this->parsimonyWarp   = parsimonyWarp   ; }
  void setEsprStopProp ( double esprStopProp ) {this->esprStopProp = esprStopProp ; }

  void validateRunParams(); 	// TODO  

private: 
  void initTrees(vector<TreeAln*> &trees, const CommandLine &cl ); 
  bool convergenceDiagnostic(); 

  vector<CoupledChains> runs; 


  // move options
  double esprStopProp; 
  double parsimonyWarp;   
  int guidedRadius; 

  ParallelSetup pl; 

  double initTime; 

  BlockRunParameters runParams;
};  


#endif
