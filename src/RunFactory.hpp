/**
   @file RunFactory

   @brief an essential class that uses all available information to
   create the exact setup of the run in terms of priors, proposals and
   so on


*/ 


#ifndef _PROPOSALFACTORY
#define _PROPOSALFACTORY

#include <vector>

#include "ProposalSet.hpp"
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
  /** 
      @brief configures the proposals 
  */ 
  std::vector<std::unique_ptr<AbstractProposal> >  
  produceProposals(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, 
		   const BlockParams& partitionParams, const TreeAln &traln, 
		   const unique_ptr<LikelihoodEvaluator> &eval, bool componentWiseMH, std::vector<ProposalSet> &resultPropSet); 
  /** 
      @brief get a copy of the random variables to be integrated  
   */ 
  vector<unique_ptr<AbstractParameter> > getRandomVariables() const; 
  /** 
      @brief adds secondary parameters to proposals, if necessary (currently only branch lengths)
   */ 
  void addSecondaryParameters(AbstractProposal* proposal, const std::vector<unique_ptr<AbstractParameter> > &allParameters); 

private: 			// METHODS 
  void addStandardParameters(std::vector<std::unique_ptr<AbstractParameter> > &vars, const TreeAln &traln ); 
  void addStandardPrior(AbstractParameter* var, const TreeAln& traln ); 
  void addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<unique_ptr<AbstractParameter> > &variables); 

private: 			// ATTRIBUTES 
  vector<unique_ptr<AbstractParameter> > randomVariables; 
}; 

#endif
