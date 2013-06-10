#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"


/**
   @brief yields a set of proposls for integrating a category  
 */
void ProposalRegistry::getProposals(category_t cat, const BlockProposalConfig &config, vector<unique_ptr<AbstractProposal> >& result) const 
{
  // vector<unique_ptr<AbstractProposal> > result; 

  vector<aaMatrix_t> someMatrices; 

  // TODO: this is not terrible as a design ... but reflection of course would be cooler  
  for(int i = 0; i < NUM_PROPOSALS; ++i)
    {      
      AbstractProposal *proposal = nullptr; 
  
      proposal_type curProp = proposal_type(i); 

      switch(curProp)
	{	      
	case ST_NNI: 
	  proposal = new StatNNI(INIT_NNI_MULT);
	  break; 
	case BRANCH_LENGTHS_MULTIPLIER:	      
	  proposal = new BranchLengthMultiplier( INIT_BL_MULT) ; 
	  break; 
	case NODE_SLIDER:
	  proposal = new NodeSlider(INIT_NODE_SLIDER_MULT); 
	  break; 
	case REVMAT_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal, RevMatParameter>( INIT_RATE_SLID_WIN, "revMatSlider"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case FREQUENCY_SLIDER:
	  proposal = new PartitionProposal<SlidingProposal, FrequencyParameter>(  INIT_FREQ_SLID_WIN, "freqSlider"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 		  
	case TL_MULT:
	  proposal = new TreeLengthMultiplier(  INIT_TL_MULTI); 
	  break; 
	case E_TBR: 
	  proposal = new ExtendedTBR(  config.getEsprStopProp(), INIT_ESPR_MULT); 
	  break; 
	case E_SPR: 
	  proposal = new ExtendedSPR(  config.getEsprStopProp(), INIT_ESPR_MULT); 
	  break; 
	case PARSIMONY_SPR:	
	  proposal = new ParsimonySPR(  config.getParsimonyWarp(), INIT_ESPR_MULT); 
	  break; 
	case RATE_HET_MULTI: 
	  proposal = new PartitionProposal<MultiplierProposal,RateHetParameter>(  INIT_GAMMA_MULTI, "rateHetMulti"); 
	  proposal->setRelativeWeight(1); 
	  break; 
	case RATE_HET_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal,RateHetParameter>(  INIT_GAMMA_SLID_WIN, "rateHetSlider"); 
	  proposal->setRelativeWeight(0); 
	  break; 
	case FREQUENCY_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,FrequencyParameter>(  INIT_DIRICHLET_ALPHA, "freqDirich"); 
	  proposal->setRelativeWeight(0.5); 
	  break; 
	case REVMAT_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,RevMatParameter>( INIT_DIRICHLET_ALPHA, "revMatDirich"); 	      
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
	result.push_back(unique_ptr<AbstractProposal>(proposal)); 
    }
} 

