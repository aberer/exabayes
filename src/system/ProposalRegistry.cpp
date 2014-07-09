#include <memory>
#include <limits>
#include <unordered_map>

#include "ProposalRegistry.hpp"

#include "WeibullProposer.hpp"
#include "proposals/DivTimeProposal.hpp"
#include "proposals/LikelihoodSPR.hpp"
#include "proposals/DistributionBranchLength.hpp"
#include "proposals/ExtendedTBR.hpp" 
#include "proposals/ExtendedSPR.hpp"
#include "proposals/ParsimonySPR.hpp"
#include "proposals/StatNNI.hpp"
#include "proposals/BranchLengthMultiplier.hpp"
#include "proposals/AminoModelJump.hpp"
#include "proposals/NodeSlider.hpp"
#include "proposals/TreeLengthMultiplier.hpp"
#include "proposals/AbstractProposal.hpp"
// #include "proposals/GammaDistributionSlider.hpp"

#include "proposers/AbstractProposer.hpp"
#include "proposers/SlidingProposer.hpp"
#include "proposers/MultiplierProposer.hpp"
#include "proposers/DirichletProposer.hpp"
#include "proposers/RateDirichletProposer.hpp"
#include "proposers/RateSlidingProposer.hpp"

#include "system/extensions.hpp"

#include "proposals/ProposalType.hpp"
#include "BoundsChecker.hpp"

const double ProposalRegistry::initBranchLengthMultiplier = 1.386294; 

const double ProposalRegistry::initRateSlidingWindow = 0.15 ;

const double ProposalRegistry::initFrequencySlidingWindow = 0.2 ; 
const double ProposalRegistry::initGammaSlidingWindow = 0.75; 
const double ProposalRegistry::initSecondaryBranchLengthMultiplier = 0.098; 
const double ProposalRegistry::initTreeLengthMultiplier = 1.386294; 

const double ProposalRegistry::initDirichletAlpha = 100 ; 
const double ProposalRegistry::initDirichletProtAlpha = 50 ; 

const double ProposalRegistry::initGammaMultiplier = 0.811 ; 
const double ProposalRegistry::initNodeSliderMultiplier = 0.191 ; 


// const int ProposalRegistry::likeSprMinRadius = 1; 
// const int ProposalRegistry::likeSprMaxRadius = 3; 
// const double ProposalRegistry::likeSprWarp = 1.; 


/**
   @brief yields a set of proposls for integrating a category  
 */
std::vector<std::unique_ptr<AbstractProposal> >
ProposalRegistry::getSingleParameterProposals(Category cat, const BlockProposalConfig &config, const TreeAln &traln) const 
{
  auto result = std::vector<unique_ptr<AbstractProposal> >{} ; 

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
	  proposal = make_unique<StatNNI>(initSecondaryBranchLengthMultiplier) ;
	  break; 
	case ProposalType::BRANCH_LENGTHS_MULTIPLIER:	      
	  proposal = make_unique< BranchLengthMultiplier>( initBranchLengthMultiplier) ; 
	  break; 
	case ProposalType::NODE_SLIDER:
	  proposal = make_unique<NodeSlider> (initNodeSliderMultiplier); 
	  break; 
	case ProposalType::TL_MULT:
	  proposal = make_unique< TreeLengthMultiplier> (  ProposalRegistry::initTreeLengthMultiplier); 
	  break; 
	case ProposalType::E_TBR: 
	  {
	    auto etbrStop = config.getEtbrStopProb(); 
	    proposal = make_unique< ExtendedTBR> (  etbrStop, initSecondaryBranchLengthMultiplier); 
	  }
	  break; 
	case ProposalType::E_SPR: 
	  proposal = make_unique<ExtendedSPR> (  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case ProposalType::PARSIMONY_SPR:	
	  {
	    // decide upon radius 
	    int radius = config.getParsSPRRadius(); 
	    nat numTax = traln.getNumberOfTaxa(); 
	    if(radius == -1)
	      radius = std::floor( std::log(numTax) * 2  );
	    proposal = make_unique<ParsimonySPR> (  config.getParsimonyWarp(), initSecondaryBranchLengthMultiplier, radius ); 
	  }
	  break; 
	case ProposalType::REVMAT_SLIDER: 
	  proposal = make_unique<ParameterProposal>  (Category::SUBSTITUTION_RATES, "revMatSlider", true, 
						      make_unique<SlidingProposer>(BoundsChecker::rateMin, BoundsChecker::rateMax, true),
						      initRateSlidingWindow,1,   1e-5,1 ) ; 
	  break; 
	case ProposalType::FREQUENCY_SLIDER:
	  proposal = make_unique<ParameterProposal>  (Category::FREQUENCIES, "freqSlider", true, 
						      std::unique_ptr<SlidingProposer>( new SlidingProposer(BoundsChecker::freqMin, std::numeric_limits<double>::max(), false)), 
						      initFrequencySlidingWindow,0.5,     1e-5, 1 ) ; 
	  break; 		  
	case ProposalType::RATE_HET_MULTI: 
	  proposal = make_unique<ParameterProposal>  (Category::RATE_HETEROGENEITY, "rateHetMulti", false, 
						      std::unique_ptr<MultiplierProposer>(new MultiplierProposer(BoundsChecker::alphaMin, BoundsChecker::alphaMax)),
						      initGammaMultiplier,1,  1e-4, 1e2 ) ; 
	  break; 
	case ProposalType::RATE_HET_SLIDER: 
	  proposal = make_unique<ParameterProposal>  (Category::RATE_HETEROGENEITY, "rateHetSlider", false, 
						      std::unique_ptr<SlidingProposer>(new SlidingProposer(BoundsChecker::alphaMin, BoundsChecker::alphaMax, false)),   
						      initGammaSlidingWindow,0, 1e-4, 1e2 ) ; 
	  break; 
	case ProposalType::FREQUENCY_DIRICHLET: 
	  proposal = make_unique<ParameterProposal>  (Category::FREQUENCIES, "freqDirich", true, 
						      std::unique_ptr<DirichletProposer	>(new DirichletProposer	(BoundsChecker::freqMin, std::numeric_limits<double>::max(), false)), 
						      initDirichletAlpha,0.5, 1e-3, 1e4 ) ; 
	  break; 
	case ProposalType::REVMAT_DIRICHLET: 
	  proposal = make_unique<ParameterProposal>  (Category::SUBSTITUTION_RATES, "revMatDirich", true, 
						      std::unique_ptr<DirichletProposer	>(new DirichletProposer	 (BoundsChecker::rateMin, BoundsChecker::rateMax, true)), 
						      initDirichletAlpha,1, 1e-3, 1e4) ; 
	  break; 
	case ProposalType::LIKE_SPR: 
	  {
	    int rad =  config.getLikeSprMaxRadius(); 
	    if(rad == -1 )
	      rad = std::floor( std::log( traln.getNumberOfTaxa() ) * 2 ); 

	    proposal = make_unique<LikelihoodSPR>(rad , config.getLikeSprWarp() );
	  }
	  break; 
	case ProposalType::DIRICH_REVMAT_PER_RATE:
	  {
	    proposal = make_unique<ParameterProposal>(Category::SUBSTITUTION_RATES, "revMatDirichRate", true,
						      make_unique<RateDirichletProposer>( BoundsChecker::rateMin, BoundsChecker::rateMax),
						      initDirichletProtAlpha, 
						      // 4, 
						      1, 
						      1e-3, 1e4
						      );
	    proposal->setForProteinsOnly(true); 
	  }
	  break; 
	case ProposalType::SLIDING_REVMAT_PER_RATE:
	  {
	    proposal = make_unique<ParameterProposal>(Category::SUBSTITUTION_RATES, "revMatSliderRate", true,
						      make_unique<RateSlidingProposer>( BoundsChecker::rateMin, BoundsChecker::rateMax),
						      initRateSlidingWindow, 0.5, 0,0 // TODO !!! 
						      );
	    
	    // TODO not ready yet   
	    assert(0); 
	  }
	  break; 
	case ProposalType::AMINO_MODEL_JUMP: 
	  proposal = make_unique<AminoModelJump>();
	  break; 
	case ProposalType::BRANCH_DIST_GAMMA: 
	  proposal = make_unique<DistributionBranchLength<GammaProposer> >();
	  break;
	case ProposalType::BRANCH_SLIDER: 
	  continue; 		// TODO implement  
	  break; 
	case ProposalType::BL_DIST_WEIBULL: 
	  proposal = make_unique< DistributionBranchLength<WeibullProposer> >();
	  break; 
	case ProposalType::DIV_TIME_DIRICH:
	  proposal = make_unique<DivTimeProposal>(); 
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
