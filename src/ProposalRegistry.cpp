#include <memory>
#include <limits>
#include <unordered_map>

#include "ProposalRegistry.hpp"

#include "proposals/LikelihoodSPR.hpp"
#include "proposals/GibbsBranchLength.hpp"
#include "proposals/ExtendedTBR.hpp" 
#include "proposals/ExtendedSPR.hpp"
#include "proposals/ParsimonySPR.hpp"
#include "proposals/StatNNI.hpp"
#include "proposals/BranchLengthMultiplier.hpp"
#include "proposals/AminoModelJump.hpp"
#include "proposals/NodeSlider.hpp"
#include "proposals/TreeLengthMultiplier.hpp"
#include "proposals/AbstractProposal.hpp"

#include "proposers/AbstractProposer.hpp"
#include "proposers/SlidingProposal.hpp"
#include "proposers/MultiplierProposal.hpp"
#include "proposers/DirichletProposal.hpp"

#include "ProposalType.hpp"
#include "BoundsChecker.hpp"

const double ProposalRegistry::initBranchLengthMultiplier = 1.386294; 
const double ProposalRegistry::initRateSlidingWindow = 0.15 ;
const double ProposalRegistry::initFrequencySlidingWindow = 0.2 ; 
const double ProposalRegistry::initGammaSlidingWindow = 0.75; 
const double ProposalRegistry::initSecondaryBranchLengthMultiplier = 0.098; 
const double ProposalRegistry::initTreeLengthMultiplier = 1.386294; 
const double ProposalRegistry::initDirichletAlpha = 100 ; 
const double ProposalRegistry::initGammaMultiplier = 0.811 ; 
const double ProposalRegistry::initNodeSliderMultiplier = 0.191 ; 


const int ProposalRegistry::likeSprMinRadius = 1; 
const int ProposalRegistry::likeSprMaxRadius = 4; 
const double ProposalRegistry::likeSpWarp = 1; 


/**
   @brief yields a set of proposls for integrating a category  
 */
vector<unique_ptr<AbstractProposal> >
ProposalRegistry::getSingleParameterProposals(Category cat, const BlockProposalConfig &config, const TreeAln &traln) const 
{
  std::vector<unique_ptr<AbstractProposal> > result; 

  auto&& proposals = ProposalTypeFunc::getSingleParameterProposalsForCategory(cat ) ; 
  for(auto& p : proposals)
    {     

      double userWeight = 1; 
      if(config.wasSetByUser(p))
	{
	  userWeight = config.getProposalWeight(p); 
	  if(userWeight == 0)
	    continue; 
	} 
      else if( not ProposalTypeFunc::isReadyForProductiveUse(p)  )		
	continue;       

      auto&& proposal = std::unique_ptr<AbstractProposal>{}; 

      switch(p)
	{	      
	case ProposalType::ST_NNI: 
	  proposal = unique_ptr<AbstractProposal>(new  StatNNI(initSecondaryBranchLengthMultiplier)) ;
	  break; 
	case ProposalType::BRANCH_LENGTHS_MULTIPLIER:	      
	  proposal = unique_ptr< AbstractProposal>(new BranchLengthMultiplier( initBranchLengthMultiplier)) ; 
	  break; 
	case ProposalType::NODE_SLIDER:
	  proposal = unique_ptr<NodeSlider>( new NodeSlider(initNodeSliderMultiplier)); 
	  break; 
	case ProposalType::TL_MULT:
	  proposal = unique_ptr< TreeLengthMultiplier>( new TreeLengthMultiplier(  ProposalRegistry::initTreeLengthMultiplier)); 
	  break; 
	case ProposalType::E_TBR: 
	  proposal = unique_ptr< ExtendedTBR>( new ExtendedTBR(  config.getEtbrStopProb(), initSecondaryBranchLengthMultiplier)); 
	  break; 
	case ProposalType::E_SPR: 
	  proposal = unique_ptr<ExtendedSPR>( new ExtendedSPR(  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier)); 
	  break; 
	case ProposalType::PARSIMONY_SPR:	
	  proposal = unique_ptr<ParsimonySPR>( new ParsimonySPR(  config.getParsimonyWarp(), initSecondaryBranchLengthMultiplier)); 
	  break; 
	case ProposalType::REVMAT_SLIDER: 
	  proposal = 
	    std::unique_ptr<ParameterProposal> ( new ParameterProposal(Category::SUBSTITUTION_RATES, "revMatSlider", true, 
								       std::unique_ptr<SlidingProposal>(new SlidingProposal(BoundsChecker::rateMin, BoundsChecker::rateMax, true)),
								       initRateSlidingWindow,0.5 )) ; 
	  break; 
	case ProposalType::FREQUENCY_SLIDER:
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::FREQUENCIES, "freqSlider", true, 
									   std::unique_ptr<SlidingProposal>( new SlidingProposal(BoundsChecker::freqMin, std::numeric_limits<double>::max(), false)), 
									   initFrequencySlidingWindow,0.5 )) ; 
	  break; 		  
	case ProposalType::RATE_HET_MULTI: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::RATE_HETEROGENEITY, "rateHetMulti", false, 
									   std::unique_ptr<MultiplierProposal>(new MultiplierProposal(BoundsChecker::alphaMin, BoundsChecker::alphaMax)),
									   initGammaMultiplier,1. )) ; 
	  break; 
	case ProposalType::RATE_HET_SLIDER: 
	  proposal = 
	    unique_ptr<ParameterProposal> ( new ParameterProposal(Category::RATE_HETEROGENEITY, "rateHetSlider", false, 
								  std::unique_ptr<SlidingProposal>(new SlidingProposal(BoundsChecker::alphaMin, BoundsChecker::alphaMax, false)),   
								  initGammaSlidingWindow,0 )) ; 
	  break; 
	case ProposalType::FREQUENCY_DIRICHLET: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::FREQUENCIES, "freqDirich", true, 
									   std::unique_ptr<DirichletProposal>(new DirichletProposal(BoundsChecker::freqMin, std::numeric_limits<double>::max(), false)), 
									   initDirichletAlpha,0.5 )) ; 
	  break; 
	case ProposalType::REVMAT_DIRICHLET: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::SUBSTITUTION_RATES, "revMatDirich", true, 
									   std::unique_ptr<DirichletProposal>(new DirichletProposal (BoundsChecker::rateMin, BoundsChecker::rateMax, true)), 
									   initDirichletAlpha,0.5)) ; 
	  break; 
	case ProposalType::LIKE_SPR: 
	  proposal = unique_ptr<LikelihoodSPR>(new LikelihoodSPR(  likeSprMinRadius, likeSprMaxRadius, likeSpWarp));
	  break; 
	case ProposalType::AMINO_MODEL_JUMP: 
	  proposal = unique_ptr<AminoModelJump>( new AminoModelJump());
	  break; 
	case ProposalType::BRANCH_GIBBS: 
	  proposal = unique_ptr<GibbsBranchLength>( new GibbsBranchLength());
	  break;
	case ProposalType::BRANCH_SLIDER: 
	  continue; 		// TODO implement  
	  break; 
	default : 
	  {
	    cerr << "you did not implement case " << int(p) << " in ProposalRegistry.cpp" << endl; 
	    assert(0); 
	  }
	} 

      if(config.wasSetByUser(p))
	proposal->setRelativeWeight(userWeight);       
      result.push_back((std::move(proposal))); 	
    }

  return result; 
} 
