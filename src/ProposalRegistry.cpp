#include <memory>
#include <limits>
#include <unordered_map>

#include "ProposalRegistry.hpp"

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


/**
   @brief yields a set of proposls for integrating a category  
 */
vector<unique_ptr<AbstractProposal> >
ProposalRegistry::getSingleParameterProposals(Category cat, const BlockProposalConfig &config, const TreeAln &traln, const unique_ptr<LikelihoodEvaluator> &eval) const 
{
  vector<unique_ptr<AbstractProposal> > result; 
  
  vector<aaMatrix_t> someMatrices; // TODO 

  auto proposals = ProposalTypeFunc::getSingleParameterProposalsForCategory(cat ) ; 
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
								       std::unique_ptr<SlidingProposal>(new SlidingProposal(BoundsChecker::rateMin, BoundsChecker::rateMax)),
								       initRateSlidingWindow )) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::FREQUENCY_SLIDER:
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::FREQUENCIES, "freqSlider", true, 
									   std::unique_ptr<SlidingProposal>( new SlidingProposal(BoundsChecker::freqMin, std::numeric_limits<double>::max())), 
									   initFrequencySlidingWindow )) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 		  
	case ProposalType::RATE_HET_MULTI: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::RATE_HETEROGENEITY, "rateHetMulti", false, 
									   std::unique_ptr<MultiplierProposal>(new MultiplierProposal(BoundsChecker::alphaMin, BoundsChecker::alphaMax)),
									   initGammaMultiplier )) ; 
	  proposal->setRelativeWeight(1); 
	  break; 
	case ProposalType::RATE_HET_SLIDER: 
	  proposal = 
	    unique_ptr<ParameterProposal> ( new ParameterProposal(Category::RATE_HETEROGENEITY, "rateHetSlider", false, 
								  std::unique_ptr<SlidingProposal>(new SlidingProposal(BoundsChecker::alphaMin, BoundsChecker::alphaMax)),   
								  initGammaSlidingWindow )) ; 
	  proposal->setRelativeWeight(0); 
	  break; 
	case ProposalType::FREQUENCY_DIRICHLET: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::FREQUENCIES, "freqDirich", true, 
									   std::unique_ptr<DirichletProposal>(new DirichletProposal(BoundsChecker::freqMin, std::numeric_limits<double>::max())), 
									   initDirichletAlpha )) ; 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::REVMAT_DIRICHLET: 
	  proposal = unique_ptr<ParameterProposal> ( new ParameterProposal(Category::SUBSTITUTION_RATES, "revMatDirich", true, 
									   std::unique_ptr<DirichletProposal>(new DirichletProposal (BoundsChecker::rateMin, BoundsChecker::rateMax)), 
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



#define MULTI_PARTITION_WEIGHT 1 

// TODO remodel this entire function for productive use 
vector<unique_ptr<AbstractProposal> >  
ProposalRegistry::getMultiParameterProposals(std::vector<AbstractParameter*> params, const BlockProposalConfig &config, 
					     const TreeAln &traln, const unique_ptr<LikelihoodEvaluator> &eval)
{
  std::vector<std::unique_ptr<AbstractProposal> >  result; 

  // TODO this function should become more intelligent. Here I'll just
  // hard code that we could have proposal per category 

  std::unordered_map<Category, std::vector<AbstractParameter*>, CategoryHash > paramsPerCategory; 
  for(auto &p : params)
    {
      auto &v =  paramsPerCategory[p->getCategory()]; 
      v.push_back(p); 
    }

  for(auto &elem : paramsPerCategory)
    {
      tout << "for category " << elem.first << " we have " << std::endl; 
      for(auto &e : elem.second)
	tout << e << std::endl;       
    }
  
  // initialize a revmat proposal 
  std::unique_ptr<DirichletProposal> p(new DirichletProposal(BoundsChecker::rateMin, BoundsChecker::rateMax));

  auto proposal =  new AlignmentProposal(Category::SUBSTITUTION_RATES, 
					 "revMatDirichAll", 
					 initDirichletAlpha, 
					 paramsPerCategory[Category::SUBSTITUTION_RATES].size(), // 
					 p.get() 
					 // ,eval->clone()
); 

  // BAD =/ 
  double userWeight = 1; 
  ProposalType myType = ProposalType::DIRICH_REVMAT_ALL; 
  if(config.wasSetByUser(myType))	
    {
      userWeight = config.getProposalWeight(myType); 
      proposal->setRelativeWeight(userWeight); 
      if(userWeight == 0)
	return {}; 			// in lieu of continue 
    }
  else if( not ProposalTypeFunc::isReadyForProductiveUse(myType) )
    return {}; 


  for(auto &p : paramsPerCategory[Category::SUBSTITUTION_RATES])
    proposal->addPrimaryParameter(std::unique_ptr<AbstractParameter>(p->clone())); 
  result.emplace_back(proposal);

  return result; 
}
