#include "RunFactory.hpp"
#include "allProposals.hpp"
#include "output.h"

#include "ProposalFunctions.hpp"
#include "Parameters.hpp"

// for dummy 
#include "BlockProposalConfig.hpp"

#define TODO   0


#if 0 
void RunFactory::setupProposals(vector<Category> &proposalCategories, vector<double> proposalWeights, const PriorBelief &prior)
{

  BlockProposalConfig propConfig; 

  vector<AbstractProposal*> prop; 
  
  vector<aaMatrix_t> someMatrices; // TODO, probably extract from prior belief 

  // initialize proposals 
  for(int i = 0; i < NUM_PROPOSALS ; ++i)
    { 
      double weight = proposalWeights[proposal_type(i)]; 
      if( weight != 0)
	{
	  AbstractProposal *proposal = NULL; 
	  
	  switch(proposal_type(i))
	    {	      
	    case BRANCH_LENGTHS_MULTIPLIER:	      
	      proposal = new BranchLengthMultiplier( INIT_BL_MULT) ; 
	      break; 
	    case NODE_SLIDER:
	      proposal = new NodeSlider(  INIT_NODE_SLIDER_MULT); 
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
	      proposal = new ExtendedTBR(  propConfig.getEsprStopProp(), INIT_ESPR_MULT); 
	      break; 
	    case E_SPR: 
	      proposal = new ExtendedSPR(  propConfig.getEsprStopProp(), INIT_ESPR_MULT); 
	      break; 
	    case PARSIMONY_SPR:	
	      proposal = new ParsimonySPR(  propConfig.getParsimonyWarp(), INIT_ESPR_MULT); 
	      break; 
	    case ST_NNI: 
	      proposal = new StatNNI(  INIT_NNI_MULT); 
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
	      proposal = new RadiusMlSPR(  propConfig.getGuidedRadius() ); 
	      break; 
	    case BRANCH_COLLAPSER:
	      proposal = new BranchCollapser(); 
	      break; 
	    case AMINO_MODEL_JUMP: 
	      proposal = new AminoModelJump( someMatrices);
	      break; 
	    default : 
	      assert(0); 
	    }

	  if(not prior.categoryIsFixed(proposal->getCategory()))
	    prop.push_back(proposal);
	  else 
	    tout << "Notice: Relative weight " << weight << " was specified for proposal "   << proposal->getName()  << ". Since the prior for this category was set to fixed, this proposal will be ignored." << endl; 	      
	}
    }



#if TODO 
  // get total sum 
  double sum = 0; 
  for(auto p : prop)
    sum += p->getRelativeProbability();


  // create categories 
  vector<string> allNames = {"Topology", "BranchLengths", "Frequencies", "RevMatrix", "RateHet" }; 
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    {
      // fish out the correct proposals 
      vector<AbstractProposal*> pr; 
      double catSum = 0; 
      for(auto p : prop)
	{	  
	  if(p->getCategory() == i)
	    {
	      pr.push_back(p); 
	      catSum += p->getRelativeProbability();
	    }
	}

      if(pr.size() > 0)
	proposalCategories.push_back( Category(allNames[i], category_t(i), catSum / sum, pr )); 
    }    


#endif
  if ( isOutputProcess() )
    {
      // print some info 
      tout << "using the following moves: " << endl; 
      for(nat i = 0; i < proposalCategories.size(); ++i)
	{
	  // TODO also to info file 
      
	  tout << proposalCategories[i].getName() << " " << fixed << setprecision(2) << proposalCategories[i].getCatFreq() * 100 << "%\t" ; 
	  auto cat = proposalCategories[i]; 
	  auto p1 =  cat.getProposals(); 
#if TODO 
	  for(auto p : p1)
	    tout << "\t" << p->getName() << "(" << fixed << setprecision(2) <<  p->getRelativeProbability() * 100 << "%)" ;
	  tout << endl; 
#endif
	}
    }
}

#endif



void RunFactory::addStandardParameters(vector<RandomVariable> &vars, const TreeAln &traln )
{
  vector<bool> categoryIsActive(false, NUM_PROP_CATS);
  for(auto v : vars )
    categoryIsActive[v.getCategory()] = true; 

  nat highestId = vars[vars.size()-1].getId(); 

  // add standard stuff, if not defined yet
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    {
      category_t catIter = category_t(i) ; 
      switch(catIter)
	{	  
	case TOPOLOGY: 
	  {
	    // deafult: everything linked 
	    if(not categoryIsActive[catIter])
	      {
		RandomVariable r(catIter,highestId);
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  r.addPartition(j); 
		vars.push_back(r); 
	      }
	    break; 
	  }

	case BRANCH_LENGTHS: 
	  {
	    // default: everything is linked
	    if(not categoryIsActive[catIter])
	      {
		RandomVariable r(catIter,highestId);
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  r.addPartition(j); 
		vars.push_back(r); 
	      }
	    break; 
	  }

	case FREQUENCIES: 
	  {
	    break; 
	  }

	case SUBSTITUTION_RATES: 
	  {
	    if(not categoryIsActive[catIter] )
	      {
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  {
		    RandomVariable r(catIter, highestId);
		    r.addPartition(j); 
		    vars.push_back(r);	   
		  }
	      }
	    break;  
	  }

	case RATE_HETEROGENEITY: 
	  {
	    if(not categoryIsActive[catIter])
	      {
		for(int j = 0; j < traln.getNumberOfPartitions();++ j)
		  {
		    RandomVariable r(catIter, highestId);
		    r.addPartition(j);
		    vars.push_back(r);
		  }
	      }
	    break; 
	  }
	  
	case AA_MODEL: 
	  {
	    if( not categoryIsActive[catIter])
	      {
		for(int j = 0; j < traln.getNumberOfPartitions(); ++j)
		  {
		    pInfo* partition = traln.getPartition(i);
		    if(partition->dataType == AA_DATA)
		      {
			RandomVariable r(catIter, highestId);
			r.addPartition(j); 
			vars.push_back(r);
		      }
		  }
	      }
	    break; 
	  }
	default: assert(0); 
	}
    }
}


void RunFactory::addStandardPrior(RandomVariable &var, const TreeAln& traln )
{
  switch(var.getCategory())			// TODO such switches should be part of an object
    {
    case TOPOLOGY:  
      var.setPrior(shared_ptr<AbstractPrior>(new UniformPrior(0,0))); // TODO : proper topology prior? 
      break; 
    case BRANCH_LENGTHS: 
      var.setPrior(shared_ptr<AbstractPrior>(new ExponentialPrior(10.0)));
      break; 
    case FREQUENCIES: 
      {
	pInfo *partition = traln.getPartition(var.getPartitions()[0]);
	assert(partition->dataType == DNA_DATA || partition->dataType == AA_DATA); 

	vector<double>badHardcoded; 
	for(int i = 0; i < partition->states; ++i)
	  badHardcoded.push_back(1.); 
	var.setPrior(shared_ptr<AbstractPrior>(new DirichletPrior(badHardcoded ))); 
      }
      break; 
    case SUBSTITUTION_RATES: 
      {
	pInfo *partition = traln.getPartition(var.getPartitions()[0]);
	assert(partition->dataType == DNA_DATA); 
	
	vector<double> subAlpha = {1,1,1,1,1,1}; 
	var.setPrior(shared_ptr<AbstractPrior>(new DirichletPrior( subAlpha ))); 
      }
      break; 
    case RATE_HETEROGENEITY: 
      var.setPrior(shared_ptr<AbstractPrior>(new UniformPrior(1e-6, 200)));     
      break; 
    case AA_MODEL : 
      assert(NOT_IMPLEMENTED); 
      break; 
    default: assert(0); 
    }
}


void RunFactory::addPriorsToVariables(const TreeAln &traln,  const BlockPrior &priorInfo, vector<RandomVariable> &variables)
{
  vector<shared_ptr<AbstractPrior> > generalPriors = priorInfo.getGeneralPriors();
  vector< map<nat, shared_ptr<AbstractPrior > > >  specificPriors = priorInfo.getSpecificPriors();

  for(RandomVariable &v : variables)
    {
      auto partitionIds = v.getPartitions(); 

      map<nat, shared_ptr<AbstractPrior > >  idMap = specificPriors[ v.getCategory() ]; 
      AbstractPrior* thePrior = nullptr; 
      for(nat partId : partitionIds)	
	{	  	  
	  if(idMap.find(partId) != idMap.end()) // found 
	    {
	      if(thePrior != nullptr)
		{
		  cerr << "while setting prior for random variable " << v << ": it seems you have defined a prior for more than one partition of a set of partitions that is linked for a specific variable to estimate. Please only specify exactly one prior per variable for a set of linked partitions."  << endl; 
		  exit(1); 
		}
	      else 
		thePrior = idMap.at(partId).get(); 
	    }	 	    
	}
      
      // use a general prior, if we have not found anything
      if(thePrior == nullptr)
	thePrior = generalPriors.at(v.getCategory()).get(); 
      
      if(thePrior != nullptr)
	v.setPrior(shared_ptr<AbstractPrior>(thePrior));
      else 
	addStandardPrior(v, traln);
    }
}


void RunFactory::configureRuns(const BlockProposalConfig &propConfig, const BlockPrior &priorInfo, const BlockParams& partitionParams, const TreeAln &traln)
{
  auto vars = partitionParams.getParameters();
  addStandardParameters(vars, traln);
  addPriorsToVariables(traln, priorInfo, vars);

  cout << "Random variables: " << endl; 
  for(auto v : vars)
    cout <<  v << endl; 
  cout << endl;
}

