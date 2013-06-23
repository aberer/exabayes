/**
   @file RunFactory

   @brief an essential class that uses all available information to
   create the exact setup of the run in terms of priors, proposals and
   so on


*/ 


#ifndef _PROPOSALFACTORY
#define _PROPOSALFACTORY

#include <vector>

#include "config/BlockProposalConfig.hpp"
#include "config/BlockPrior.hpp"
#include "config/BlockParams.hpp"
#include "AbstractProposal.hpp"
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"
#include "RandomVariable.hpp"


class RunFactory
{
public: 
  RunFactory(){}
 
  void configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, vector<unique_ptr<AbstractProposal> > &proposalResult);   
  vector<RandomVariablePtr> getRandomVariables() const {return randomVariables; }

private: 
  vector<RandomVariablePtr> randomVariables; 

  void addStandardParameters(vector<RandomVariablePtr> &vars, const TreeAln &traln ); 
  void addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<RandomVariablePtr> &variables); 
  void addStandardPrior(RandomVariablePtr &var, const TreeAln& traln ); 
}; 

#endif
