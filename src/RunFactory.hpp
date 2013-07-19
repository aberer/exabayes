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
#include "parameters/AbstractParameter.hpp"


class RunFactory
{
public:  
  void configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln, vector<unique_ptr<AbstractProposal> > &proposals, shared_ptr<LikelihoodEvaluator> eval);
  vector<unique_ptr<AbstractParameter> > getRandomVariables() const
  {
    vector<unique_ptr<AbstractParameter> > result; 
    for(auto &r : randomVariables)
      result.push_back(std::unique_ptr<AbstractParameter>(r->clone()));
    return result; 
  }

private: 
  vector<unique_ptr<AbstractParameter> > randomVariables; 

  void addStandardParameters(vector<unique_ptr<AbstractParameter> > &vars, const TreeAln &traln ); 
  void addStandardPrior(AbstractParameter* var, const TreeAln& traln ); 
  void addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<unique_ptr<AbstractParameter> > &variables); 
}; 

#endif
