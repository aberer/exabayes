/**
   @file RunFactory

   @brief an essential class that uses all available information to
   create the exact setup of the run in terms of priors, proposals and
   so on


*/ 


#ifndef _PROPOSALFACTORY
#define _PROPOSALFACTORY

#include <vector>

#include "proposals/ProposalSet.hpp"
#include "config/BlockProposalConfig.hpp"
#include "config/BlockPrior.hpp"
#include "config/BlockParams.hpp"
#include "proposals/AbstractProposal.hpp"
#include "model/TreeAln.hpp"
#include "GlobalVariables.hpp"
#include "parameters/AbstractParameter.hpp"


class RunFactory
{
public:  
  /** 
      @brief configures the proposals 
  */ 
  std::tuple<std::vector<std::unique_ptr<AbstractProposal> >, std::vector<ProposalSet> >
  produceProposals(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, 
		   std::vector<std::unique_ptr<AbstractParameter>>  & params, 
		   const TreeAln &traln, bool componentWiseMH  ); 
  /** 
      @brief get a copy of the random variables to be integrated  
   */ 
  std::vector<std::unique_ptr<AbstractParameter> > getParametersToIntegrate() const; 
  void addStandardParameters(std::vector<std::unique_ptr<AbstractParameter> > &vars, const TreeAln &traln ) const; 
private: 			// METHODS 
  void addStandardPrior(AbstractParameter* var, const TreeAln& traln ); 
  void addPriorsToParameters(const TreeAln &traln,  const BlockPrior &priorInfo, vector<unique_ptr<AbstractParameter> > &variables); 
  /** 
      @brief adds secondary parameters to proposals, if necessary (currently only branch lengths)
  */ 
  void addSecondaryParameters(AbstractProposal* proposal, const std::vector<unique_ptr<AbstractParameter> > &allParameters); 
}; 

#endif
