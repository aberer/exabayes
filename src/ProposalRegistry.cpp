#include "ProposalRegistry.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"


/**
   @brief yields a set of proposls for integrating a category  
 */
vector<shared_ptr<AbstractProposal> > ProposalRegistry::getProposals(category_t cat, const BlockProposalConfig &config) const 
{
  vector<shared_ptr<AbstractProposal> > result; 

  vector<aaMatrix_t> someMatrices; 

  // TODO: this is not terrible as a design ... but reflection of course would be cooler  
  for(int i = 0; i < NUM_PROPOSALS; ++i)
    {      
      AbstractProposal *proposal = nullptr; 
  
      switch(proposal_type(i))
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
	  break; 
	case FREQUENCY_SLIDER:
	  proposal = new PartitionProposal<SlidingProposal, FrequencyParameter>(  INIT_FREQ_SLID_WIN, "freqSlider"); 
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
	  break; 
	case RATE_HET_SLIDER: 
	  proposal = new PartitionProposal<SlidingProposal,RateHetParameter>(  INIT_GAMMA_SLID_WIN, "rateHetSlider"); 
	  break; 
	case FREQUENCY_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,FrequencyParameter>(  INIT_DIRICHLET_ALPHA, "freqDirich"); 
	  break; 
	case REVMAT_DIRICHLET: 
	  proposal = new PartitionProposal<DirichletProposal,RevMatParameter>( INIT_DIRICHLET_ALPHA, "revMatDirich"); 	      
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

      if(proposal->getCategory() != cat)
	delete proposal; 
      else 
	result.push_back(shared_ptr<AbstractProposal>(proposal)); 
    }
  return result;
} 


void ProposalRegistry::updateProposalWeights(const BlockProposalConfig &propConfig) const
{
  auto weights =  propConfig.getUserProposalWeights();

  double *dummy = nullptr; 	// specify the variable in your proposal, to enable the proposal
  

  // maps the type to the location of the relative weight variable. Just copy the scheme 
  map<proposal_type,double*> theMap =   
    {
      {   ST_NNI , &(StatNNI::relativeWeight)} ,
      {	  E_SPR, &(ExtendedSPR::relativeWeight) },
      {	  E_TBR, &(ExtendedTBR::relativeWeight) },
      {	  PARSIMONY_SPR, &(ParsimonySPR::relativeWeight) },
      {	  GUIDED_SPR, &(RadiusMlSPR::relativeWeight) },
      {	  REVMAT_SLIDER, &(PartitionProposal<SlidingProposal,RevMatParameter>::relativeWeight) },
      {	  REVMAT_DIRICHLET, &(PartitionProposal<DirichletProposal,RevMatParameter>::relativeWeight) },
      {	  RATE_HET_SLIDER, &(PartitionProposal<SlidingProposal,RateHetParameter>::relativeWeight) },
      {	  RATE_HET_MULTI, &(PartitionProposal<MultiplierProposal,RateHetParameter>::relativeWeight) },
      {	  FREQUENCY_SLIDER, &(PartitionProposal<SlidingProposal,FrequencyParameter>::relativeWeight) },
      {	  FREQUENCY_DIRICHLET, &(PartitionProposal<DirichletProposal,FrequencyParameter>::relativeWeight) },
      {	  TL_MULT, &(TreeLengthMultiplier::relativeWeight) },
      {	  BRANCH_COLLAPSER, &(BranchCollapser::relativeWeight) },
      {	  NODE_SLIDER, &(NodeSlider::relativeWeight) },
      {	  BRANCH_LENGTHS_MULTIPLIER, &(BranchLengthMultiplier::relativeWeight) }, 
      {   BRANCH_SLIDER , dummy}, 
      {   UPDATE_SINGLE_BL_GUIDED, dummy},
      {   AMINO_MODEL_JUMP, dummy}
    } ; 

  for(nat i = 0; i < NUM_PROPOSALS; ++i)
    {
      if(theMap[proposal_type(i)] != nullptr)
	{
	  *(theMap[proposal_type(i)]) = weights[i]; 
	}
    }
}
