#ifndef _PROPOSAL_REGISTRY_H
#define _PROPOSAL_REGISTRY_H

#include <vector>
#include <memory>

#include "GibbsBranchLength.hpp"
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
#include "config/BlockProposalConfig.hpp"

class ProposalRegistry
{
public: 
  void getProposals(Category cat, const BlockProposalConfig &config, 
		    vector<unique_ptr<AbstractProposal> > &result, 
		    const TreeAln &traln, shared_ptr<LikelihoodEvaluator> eval) const ; 
  
  static const double initFrequencySlidingWindow ; 
  static const  double initBranchLengthMultiplier; 
  static const  double initRateSlidingWindow; 
  static const double initGammaSlidingWindow; 
  static const double initSecondaryBranchLengthMultiplier; 
  static const double initDirichletAlpha; 
  static const double initTreeLengthMultiplier ; 
  static const double initGammaMultiplier; 
  static const double initNodeSliderMultiplier; 


private: 

}; 


#endif
