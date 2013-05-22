#ifndef _PROPOSAL_REGISTRY_H
#define _PROPOSAL_REGISTRY_H

#include <vector>
#include <memory>

#include "ExtendedTBR.hpp"
#include "ExtendedSPR.hpp"
#include "ParsimonySPR.hpp"
#include "StatNNI.hpp"
#include "RadiusMlSPR.hpp"
#include "BranchLengthMultiplier.hpp"
#include "BranchCollapser.hpp"
#include "AminoModelJump.hpp"
#include "NodeSlider.hpp"
#include "TreeLengthMultiplier.hpp"
#include "PartitionProposal.hpp"
#include "AbstractProposal.hpp"
#include "GlobalVariables.hpp"
#include "BlockProposalConfig.hpp"

class ProposalRegistry
{
public: 
  void  getProposals(category_t cat, const BlockProposalConfig &config, vector<unique_ptr<AbstractProposal> > &result) const ;
  void updateProposalWeights(const BlockProposalConfig &propConfig) const; 

private: 

}; 


#endif
