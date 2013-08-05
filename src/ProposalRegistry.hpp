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
#include "AbstractProposal.hpp"
#include "GlobalVariables.hpp"
#include "config/BlockProposalConfig.hpp"
#include "ParameterProposal.hpp"
#include "proposals/AlignmentProposal.hpp"

class ProposalRegistry
{
public: 
  /** 
      @brief get all proposals that integrate over a single parameter   
  */ 
  vector<unique_ptr<AbstractProposal> >
  getSingleParameterProposals(Category cat, const BlockProposalConfig &config, const TreeAln &traln, const unique_ptr<LikelihoodEvaluator> &eval) const ; 
  /** 
      @brief get proposals that integrate over multiple parameters 
   */ 
  vector<unique_ptr<AbstractProposal> >  
  getMultiParameterProposals(std::vector<AbstractParameter*> params, const BlockProposalConfig &config, const TreeAln &traln, const unique_ptr<LikelihoodEvaluator> &eval); 

  static const double initFrequencySlidingWindow ; 
  static const double initBranchLengthMultiplier; 
  static const double initRateSlidingWindow; 
  static const double initGammaSlidingWindow; 
  static const double initSecondaryBranchLengthMultiplier; 
  static const double initDirichletAlpha; 
  static const double initTreeLengthMultiplier ; 
  static const double initGammaMultiplier; 
  static const double initNodeSliderMultiplier; 

private: 

}; 


#endif
