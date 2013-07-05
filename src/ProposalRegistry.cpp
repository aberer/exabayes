#include <memory>
#include <limits>

#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
// #include "Parameters.hpp"
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


/**
   @brief yields a set of proposls for integrating a category  
 */
void ProposalRegistry::getProposals(Category cat, const BlockProposalConfig &config, vector<unique_ptr<AbstractProposal> > &result, const TreeAln &traln, shared_ptr<LikelihoodEvaluator> eval) const 
{
  vector<aaMatrix_t> someMatrices; 

  auto proposals =  ProposalTypeFunc::getProposalsForCategory(cat); 
  for(auto p : proposals)
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

      unique_ptr<AbstractProposal> proposal; 

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
	  proposal = unique_ptr< ExtendedTBR>( new ExtendedTBR(  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier)); 
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
								       make_shared<SlidingProposal>(BoundsChecker::rateMin, BoundsChecker::rateMax),
								       initRateSlidingWindow )) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::FREQUENCY_SLIDER:
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::FREQUENCIES, "freqSlider", true, 
									   make_shared<SlidingProposal>(BoundsChecker::freqMin, std::numeric_limits<double>::max()), 
									   initFrequencySlidingWindow )) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 		  
	case ProposalType::RATE_HET_MULTI: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::RATE_HETEROGENEITY, "rateHetMulti", false, 
									   make_shared<MultiplierProposal>(BoundsChecker::alphaMin, BoundsChecker::alphaMax),
									   initGammaMultiplier )) ; 
	  proposal->setRelativeWeight(1); 
	  break; 
	case ProposalType::RATE_HET_SLIDER: 
	  proposal = 
	    unique_ptr<ParameterProposal> ( new ParameterProposal(Category::RATE_HETEROGENEITY, "rateHetSlider", false, 
								  make_shared<SlidingProposal>(BoundsChecker::alphaMin, BoundsChecker::alphaMax),   
								  initGammaSlidingWindow )) ; 
	  proposal->setRelativeWeight(0); 
	  break; 
	case ProposalType::FREQUENCY_DIRICHLET: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::FREQUENCIES, "freqDirich", true, 
									   make_shared<DirichletProposal>(BoundsChecker::freqMin, std::numeric_limits<double>::max()), 
									   initDirichletAlpha )) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::REVMAT_DIRICHLET: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::SUBSTITUTION_RATES, "revMatDirich", true, 
									   make_shared<DirichletProposal>(BoundsChecker::rateMin, BoundsChecker::rateMax), 
									   initDirichletAlpha)) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::GUIDED_SPR:
	  proposal = unique_ptr<RadiusMlSPR>( new RadiusMlSPR(  config.getGuidedRadius() )); 
	  break; 
	case ProposalType::BRANCH_COLLAPSER:
	  proposal = unique_ptr<BranchCollapser>( new BranchCollapser()); 
	  break; 
	case ProposalType::AMINO_MODEL_JUMP: 
	  proposal = unique_ptr<AminoModelJump>( new AminoModelJump(someMatrices));
	  break; 
	case ProposalType::BRANCH_GIBBS: 
	  proposal = unique_ptr<GibbsBranchLength>( new GibbsBranchLength(eval));
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
} 

