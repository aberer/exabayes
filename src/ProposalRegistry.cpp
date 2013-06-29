#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"

// #define INIT_BL_SLID_WIN  0.075

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

  // TODO: this is not terrible as a design ... but reflection of course would be cooler  
  for(int i = 0; i < NUM_PROPOSALS; ++i)
    {
      AbstractProposal *proposal = nullptr; 
  
      proposal_type curProp = proposal_type(i); 

      switch(curProp)
	{	      
	case ST_NNI: 
	  proposal = new StatNNI(initSecondaryBranchLengthMultiplier);
	  break; 
	case BRANCH_LENGTHS_MULTIPLIER:	      
	  proposal = new BranchLengthMultiplier( initBranchLengthMultiplier) ; 
	  break; 
	case NODE_SLIDER:
	  proposal = new NodeSlider(initNodeSliderMultiplier); 
	  break; 
	case REVMAT_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal, RevMatParameter>( initRateSlidingWindow, "revMatSlider"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case FREQUENCY_SLIDER:
	  proposal = new PartitionProposal<SlidingProposal, FrequencyParameter>(  initFrequencySlidingWindow, "freqSlider"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 		  
	case TL_MULT:
	  proposal = new TreeLengthMultiplier(  ProposalRegistry::initTreeLengthMultiplier); 
	  break; 
	case E_TBR: 
	  proposal = new ExtendedTBR(  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case E_SPR: 
	  proposal = new ExtendedSPR(  config.getEsprStopProp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case PARSIMONY_SPR:	
	  proposal = new ParsimonySPR(  config.getParsimonyWarp(), initSecondaryBranchLengthMultiplier); 
	  break; 
	case RATE_HET_MULTI: 
	  proposal = new PartitionProposal<MultiplierProposal,RateHetParameter>(  initGammaMultiplier, "rateHetMulti"); 
	  proposal->setRelativeWeight(1); 
	  break; 
	case RATE_HET_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal,RateHetParameter>(  initGammaSlidingWindow, "rateHetSlider"); 
	  proposal->setRelativeWeight(0); 
	  break; 
	case FREQUENCY_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,FrequencyParameter>(  initDirichletAlpha, "freqDirich"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case REVMAT_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,RevMatParameter>( initDirichletAlpha, "revMatDirich"); 	      
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case GUIDED_SPR:
	  proposal = new RadiusMlSPR(  config.getGuidedRadius() ); 
	  break; 
	case BRANCH_COLLAPSER:
	  proposal = new BranchCollapser(); 
	  break; 
	case AMINO_MODEL_JUMP: 
	  proposal = new AminoModelJump(someMatrices);
	  break; 
	case BRANCH_GIBBS: 
	  proposal = new GibbsBranchLength(eval);
	  break;
	case UPDATE_SINGLE_BL_GUIDED: 
	case BRANCH_SLIDER: 
	  continue; 		// TODO implement  
	  break; 
	default : 
	  {
	    cerr << "you did not implement case " << i << endl; 
	    assert(0); 
	  }
	}
      
      if(config.wasSetByUser(curProp) )
	proposal->setRelativeWeight(config.getProposalWeight(curProp)); 

      if(proposal->getCategory() != cat || proposal->getRelativeWeight() == 0 )
	delete proposal; 
      else 
	result.push_back(ProposalPtr(proposal)); 
    }
} 

