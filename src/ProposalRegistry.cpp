#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"
#include "ProposalType.hpp"


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
void ProposalRegistry::getProposals(Category cat, const BlockProposalConfig &config, vector<ProposalPtr>& result, const TreeAln &traln, LikelihoodEvaluatorPtr &eval) const 
{
  vector<aaMatrix_t> someMatrices; 

  auto proposals =  ProposalTypeFunc::getProposalsForCategory(cat); 
  for(auto p : proposals)
    {     
      AbstractProposal *proposal = nullptr; 

      switch(p)
	{	      
	case ProposalType::ST_NNI: 
	  proposal = new StatNNI(initSecondaryBranchLengthMultiplier);
	  break; 
	case ProposalType::BRANCH_LENGTHS_MULTIPLIER:	      
	  proposal = new BranchLengthMultiplier( initBranchLengthMultiplier) ; 
	  break; 
	case ProposalType::NODE_SLIDER:
	  proposal = new NodeSlider(initNodeSliderMultiplier); 
	  break; 
	case ProposalType::REVMAT_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal, RevMatParameter>( initRateSlidingWindow, "revMatSlider"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::FREQUENCY_SLIDER:
	  proposal = new PartitionProposal<SlidingProposal, FrequencyParameter>(  initFrequencySlidingWindow, "freqSlider"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 		  
	case ProposalType::TL_MULT:
	  proposal = new TreeLengthMultiplier(  ProposalRegistry::initTreeLengthMultiplier); 
	  break; 
	case ProposalType::E_TBR: 
	  proposal = new ExtendedTBR(  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case ProposalType::E_SPR: 
	  proposal = new ExtendedSPR(  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case ProposalType::PARSIMONY_SPR:	
	  proposal = new ParsimonySPR(  config.getParsimonyWarp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case ProposalType::RATE_HET_MULTI: 
	  proposal = new PartitionProposal<MultiplierProposal,RateHetParameter>(  initGammaMultiplier, "rateHetMulti"); 
	  proposal->setRelativeWeight(1); 
	  break; 
	case ProposalType::RATE_HET_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal,RateHetParameter>(  initGammaSlidingWindow, "rateHetSlider"); 
	  proposal->setRelativeWeight(0); 
	  break; 
	case ProposalType::FREQUENCY_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,FrequencyParameter>(  initDirichletAlpha, "freqDirich"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::REVMAT_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,RevMatParameter>( initDirichletAlpha, "revMatDirich"); 	      
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case ProposalType::GUIDED_SPR:
	  proposal = new RadiusMlSPR(  config.getGuidedRadius() ); 
	  break; 
	case ProposalType::BRANCH_COLLAPSER:
	  proposal = new BranchCollapser(); 
	  break; 
	case ProposalType::AMINO_MODEL_JUMP: 
	  proposal = new AminoModelJump(someMatrices);
	  break; 
	case ProposalType::BRANCH_GIBBS: 
	  proposal = new GibbsBranchLength(eval);
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
	proposal->setRelativeWeight(config.getProposalWeight(p)); 

      if(proposal->getRelativeWeight() == 0 )
	delete proposal; 
      else 
	result.push_back(ProposalPtr(proposal)); 
    }
} 

