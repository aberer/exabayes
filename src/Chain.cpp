#include <sstream>

#include "Chain.hpp"
#include "Topology.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "output.h"
#include "GlobalVariables.hpp"
#include "BipartitionHash.hpp"
#include "adapters.h"
#include "proposals.h"
#include "ExtendedTBR.hpp"
#include "ExtendedSPR.hpp"
#include "WrappedProposal.hpp"
#include "ParsimonySPR.hpp" 
#include "eval.h"
#include "TreeLengthMultiplier.hpp"
#include "StatNNI.hpp"
#include "PartitionProposal.hpp"
#include "ProposalFunctions.hpp"
#include "Parameters.hpp"


// meh =/ 
extern bool isNewProposal[NUM_PROPOSALS]; 

#include <sstream>

Chain::Chain(randKey_t seed, int id, int _runid, TreeAln* _traln, string runName )  
  : traln(_traln)
  , couplingId(id)
  , currentGeneration(0)
  , hastings(1)
  , runid(_runid)
{
  chainRand = new Randomness(seed.v[0]);

  string workdir = ""; 
  assert(0); 
  
  if(id == 0)
    {
      stringstream tNameBuilder,
	pNameBuilder;
      tNameBuilder << workdir << PROGRAM_NAME << "_topologies." << runName << "."  << runid ; 
      pNameBuilder << workdir << PROGRAM_NAME << "_parameters." << runName << "."  << runid  ; 

      /* todo binary Chain file?  */
      topologyFile = fopen(tNameBuilder.str().c_str(), "w"); 
      outputParamFile = fopen(pNameBuilder.str().c_str(), "w");   
      initializeOutputFiles(this);
    }
  
  assert(0);
#if 0 
  setupProposals(initParams); 
#endif
  
  initParamDump(); 

  for(int j = 0; j < traln->getNumberOfPartitions(); ++j)
    traln->initRevMat(j);

  evaluateFullNoBackup(this);   
  
  tree *tr = traln->getTr();
  Randomness *rand = getChainRand(); 
  stringstream ss; 
  ss << *rand ;	  	  
  printInfo("init lnl=%f\tTL=%f\tseeds=%s\n", traln->getTr()->likelihood, branchLengthToReal(tr, traln->getTreeLength()), ss.str().c_str()); 

  saveTreeStateToChain(); 
}



/**
   @brief returns the inverse temperature for this chain
 */ 
double Chain::getChainHeat()
{
  double tmp = 1. + ( deltaT * couplingId ) ; 
  double inverseHeat = 1.f / tmp; 
  assert(couplingId == 0 || inverseHeat < 1.); 
  return inverseHeat; 
}


void Chain::setupProposals( initParamStruct *initParams)
{
  vector<proposalFunction*> legProp; 
  vector<AbstractProposal*> prop; 
  
  // initialize proposals 
  for(int i = 0; i < NUM_PROPOSALS ; ++i)
    { 
      if( NOT isNewProposal[i])
	{
	  proposalFunction *pf = NULL; 
	  initProposalFunction(proposal_type(i), initParams, &pf); 
	  if(pf != (proposalFunction*) NULL)
	    prop.push_back(new WrappedProposal( pf, this)); 
	}
      else 
	{
	  double weight = initParams->initWeights[proposal_type(i)]; 
	  if( weight != 0)
	    {
	      switch(proposal_type(i))
		{
		case UPDATE_MODEL: 
		  prop.push_back(new PartitionProposal<SlidingProposal, RevMatParameter>(this,weight, INIT_RATE_SLID_WIN, "revMatSlider"));
		  break; 
		case FREQUENCY_SLIDER:
		  prop.push_back(new PartitionProposal<SlidingProposal, FrequencyParameter>(this, weight, INIT_FREQ_SLID_WIN, "freqSlider"));
		  break; 		  
		case TL_MULT:
		  prop.push_back(new TreeLengthMultiplier(this, weight, INIT_TL_MULTI));
		  break; 
		case E_TBR: 
		  prop.push_back(new ExtendedTBR(this, weight, initParams->eSprStopProb, INIT_ESPR_MULT)); 
		  break; 
		case E_SPR: 
		  prop.push_back(new ExtendedSPR(this, weight, initParams->eSprStopProb, INIT_ESPR_MULT)); 
		  break; 
		case PARSIMONY_SPR:	
		  prop.push_back(new ParsimonySPR(this, weight, initParams->parsWarp, INIT_ESPR_MULT)); 
		  break; 
		case ST_NNI: 
		  prop.push_back(new StatNNI(this, weight, INIT_NNI_MULT)); 
		  break; 
		case GAMMA_MULTI: 
		  prop.push_back(new PartitionProposal<MultiplierProposal,RateHetParameter>(this, weight, INIT_GAMMA_MULTI, "rateHetMulti")); 
		  break; 
		case UPDATE_GAMMA: 
		  prop.push_back(new PartitionProposal<SlidingProposal,RateHetParameter>(this, weight, INIT_GAMMA_SLID_WIN, "rateHetSlider")); 
		  break; 
		case UPDATE_GAMMA_EXP: 
		  prop.push_back(new PartitionProposal<ExponentialProposal,RateHetParameter>(this, weight, 0, "rateHetExp")); 
		  break; 
		case UPDATE_FREQUENCIES_DIRICHLET: 
		  prop.push_back(new PartitionProposal<DirichletProposal,FrequencyParameter>(this, weight, INIT_DIRICHLET_ALPHA, "freqDirich")); 
		  break; 
		case UPDATE_MODEL_DIRICHLET: 
		  prop.push_back(new PartitionProposal<DirichletProposal,RevMatParameter>(this,weight, INIT_DIRICHLET_ALPHA, "revMatDirich"));
		  break; 
		default : 
		  assert(0); 
		}
	    }
	}  
    }

  // get total sum 
  double sum = 0; 
  for(int i = 0; i < NUM_PROPOSALS; ++i)    
    sum += initParams->initWeights[i]; 

  // create categories 
  vector<string> allNames = {"", "Topology", "BranchLengths", "Frequencies", "RevMatrix", "RateHet" }; 
  for(int i = 1; i < NUM_PROP_CATS+1; ++i)
    {
      // fish out the correct proposals 
      vector<proposalFunction*> legPr; 
      vector<AbstractProposal*> pr; 
      double catSum = 0; 
      for(auto p : legProp)
	if(p->category == i)
	  {
	    legPr.push_back(p); 
	    catSum += p->relativeWeight; 
	  }
      for(auto p : prop)
	if(p->getCategory() == i)
	  {
	    pr.push_back(p); 
	    catSum += p->getRelativeProbability();
	  }

      if(pr.size() > 0 || legPr.size() > 0)
	proposalCategories.push_back( Category(allNames[i], category_t(i), catSum / sum, pr )); 
    }    


  if ( isOutputProcess() 
       && couplingId == 0
       && runid == 0)
    {
      // print some info 
      cout << "using the following moves: " << endl; 
      for(nat i = 0; i < proposalCategories.size(); ++i)
	{
	  // TODO also to info file 
      
	  cout << proposalCategories[i].getName() << " " << fixed << setprecision(2) << proposalCategories[i].getCatFreq() * 100 << "%\t" ; 
	  auto cat = proposalCategories[i]; 
	  auto p1 =  cat.getProposals(); 
	  for(auto p : p1)
	    cout << "\t" << p->getName() << "(" << fixed << setprecision(2) <<  p->getRelativeProbability() * 100 << "%)" ;
	  cout << endl; 
	}
    }
}



/**
   @brief assigns the entire chain state from rhs to this chain

   This includes topology, parameters, all proposal parameters 
 */ 
Chain& Chain::operator=(Chain& rhs)
{
  TreeAln &thisTraln = *traln; 
  TreeAln &rhsTraln = *(rhs.traln); 

  thisTraln = rhsTraln; 

  // TODO should also copy proposals, but that is currently not used

  dump.topology->saveTopology(thisTraln);  
  return *this; 
}


void Chain::saveTreeStateToChain()
{
  dump.topology->saveTopology(*(traln));

  /* save model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = dump.infoPerPart + i; 
      pInfo *partition = traln->getPartition(i); 
      
      info->alpha =  partition->alpha; 
      
      memcpy(info->substRates, partition->substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, partition->frequencies, 4 * sizeof(double)); 
    }  
}




void Chain::applyChainStateToTree()
{
  /* TODO enable multi-branch    */
  assert(traln->getNumBranches() == 1); 

  dump.topology->restoreTopology(*(traln));

  /* restore model parameters */
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = dump.infoPerPart + i ; 

      pInfo *partition = traln->getPartition( i);
      partition->alpha = info->alpha;

      memcpy(partition->substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(partition->frequencies, info->frequencies, 4 * sizeof(double));
      
      traln->initRevMat(i);
      traln->discretizeGamma(i);	 
    } 

  evaluateFullNoBackup(this); 
}




void Chain::debug_printAccRejc(AbstractProposal *prob, bool accepted, double lnl) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  if(isOutputProcess())
    {
      if(accepted)
	printInfo( "ACC\t");   
      else 
	printInfo("rej\t");   	  
      printf("%s %g\n" ,prob->getName().c_str() , lnl); 
    }
#endif
}




AbstractProposal* Chain::drawProposalFunction()
{    
  double r = chainRand->drawRandDouble01();  
  for(auto c : proposalCategories)
    {
      double ref = c.getCatFreq(); 
      if(r <= ref )
	return c.drawProposal(*chainRand);
      else 
	r -= ref; 
    }
  
  assert(0); 
  return NULL; 
}



void Chain::step()
{
  this->currentGeneration++; 

  this->restorer->resetRestorer(*traln);

  // inform the rng that we produce random numbers for generation x  
  chainRand->rebase(currentGeneration);

  tree *tr = this->traln->getTr();   

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = getChainHeat();

  AbstractProposal *pfun = drawProposalFunction();
 
  /* reset proposal ratio  */
  this->hastings = 1; 
  
  double oldPrior = prior.getLogProb();
  pfun->applyToState(*traln, prior, hastings, *chainRand);
  pfun->evaluateProposal(*traln, prior);
  
  double priorRatio = prior.getLogProb() - oldPrior;
  double lnlRatio = tr->likelihood - prevLnl; 

  double testr = chainRand->drawRandDouble01();
  double acceptance = 
    exp(( priorRatio   + lnlRatio) * myHeat) 
    * hastings;
  
  bool wasAccepted  = testr < acceptance; 

  // printInfo("prior=%f\tlnl=%f\tacc=%f\ttestr=%f\n", priorRatio, lnlRatio, acceptance, testr); 
  
  debug_printAccRejc( pfun, wasAccepted, tr->likelihood); 


#ifdef VERIFY_LNL_SUPER_EXPENSIVE  
  // TEST
  double lnlBefore = tr->likelihood; 
  evaluateFullNoBackup(this); 
  assert( ( this->traln->getTr()->likelihood  - lnlBefore ) < ACCEPTED_LIKELIHOOD_EPS); 
  // END
#endif


  if(wasAccepted)
    {
      pfun->accept();      
      expensiveVerify(this);
    }
  else
    {
      pfun->resetState(*traln, prior);
      pfun->reject();

      // TODO maybe not here =/ 
      this->restorer->restoreArrays(*traln); // restores the previous tree state 
    }

#ifdef VERIFY_LNL_SUPER_EXPENSIVE  
  // {
  //   double lnlBefore = tr->likelihood;  
  //   evaluateFullNoBackup(this);
  //   assert( ( this->traln->getTr()->likelihood  - lnlBefore ) < ACCEPTED_LIKELIHOOD_EPS); 
  // }
#endif

  debug_checkTreeConsistency(this->traln->getTr());
 
  if(this->tuneFrequency <  pfun->getNumCallSinceTuning() )
    pfun->autotune();


#if 0 
  /* autotuning for proposal parameters. With increased parallelism
     this will become more complicated.  */
  if( globals.tuneFreq > 0 && this->currentGeneration % globals.tuneFreq == globals.tuneFreq - 1 )
    {
      for(int i = 0; i < this->numProposals; ++i)
	{
	  proposalFunction *pf = this->proposals[i]; 
	  
	  if(pf->autotune)	/* only, if we set this   */
	    {
#ifdef TUNE_ONLY_IF_ENOUGH
	      if(pf->sCtr.lAcc + pf->sCtr.lRej < globals.tuneFreq )
		continue; 
#endif
	      pf->autotune(this, pf);
	    }
	}
    }
#endif


}





void Chain::initParamDump()
{  
  tree *tr = traln->getTr();
  
  dump.topology = new Topology(tr->mxtips); 
  dump.infoPerPart = (perPartitionInfo*)exa_calloc(traln->getNumberOfPartitions(), sizeof(perPartitionInfo));
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *p =  dump.infoPerPart + i ; 
      p->alpha = 0.5; 
      for(int j = 0; j < 6; ++j)
	p->substRates[j] = 0.5; 
      p->substRates[5] = 1; 
      for(int j = 0; j < 4; ++j)
	p->frequencies[j] = 0.25;       
    }
}


/**
   @brief prints a message with associated chain/run/heat information. 

   Please always use this, when possible, it also assures that the
   message is only printed once.
 */ 
void Chain::printInfo(const char *format, ...)
{  
  if( isOutputProcess())
    {
      printf("[run %d / heat %d / gen %d] ", this->runid, this->couplingId, this->currentGeneration); 
      va_list args;
      va_start(args, format);     
      vprintf(format, args );
      va_end(args);
    }
}




void Chain::clarifyOwnership()
{
  for(auto c : proposalCategories)
    for(auto p : c.getProposals())
      p->setOwningChain(this); 
}

