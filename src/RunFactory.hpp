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
 
  void configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, vector<unique_ptr<AbstractProposal> > &proposals, shared_ptr<LikelihoodEvaluator> eval);
  vector<shared_ptr<RandomVariable> > getRandomVariables() const {return randomVariables; }

private: 
  vector<shared_ptr<RandomVariable> > randomVariables; 

  void addStandardParameters(vector<shared_ptr<RandomVariable> > &vars, const TreeAln &traln ); 
  void addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<shared_ptr<RandomVariable> > &variables); 
  void addStandardPrior(RandomVariable* var, const TreeAln& traln ); 
}; 

#endif
