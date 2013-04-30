#include <sstream>

#include "Chain.hpp"
#include "Topology.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "Randomness.hpp"
#include "output.h"
#include "globals.h"
#include "BipartitionHash.hpp"
#include "adapters.h"
#include "proposals.h"
#include "topology-utils.h"


Chain::Chain(randKey_t seed, int id, int _runid, TreeAln* _traln, initParamStruct *initParams)  
  : traln(_traln)
  , couplingId(id)
  , currentGeneration(0)
  , priorProb(1)
  , hastings(1)
  , runid(_runid)
{
  chainRand = new Randomness(seed.v[0]);
  assert(id < gAInfo.numberCoupledChains); 
  
  if(id == 0)
    {
      char tName[1024],
	pName[1024] ; 
  
      sprintf(tName, "%s%s_topologies.%s.%d", workdir, PROGRAM_NAME, run_id, runid); 
      sprintf(pName, "%s%s_parameters.%s.%d", workdir, PROGRAM_NAME, run_id, runid); 
  
      /* todo binary Chain file?  */
      topologyFile = fopen(tName, "w"); 
      outputParamFile = fopen(pName, "w");   
      initializeOutputFiles(this);
    }

  categoryWeights = (double*)exa_calloc(NUM_PROP_CATS, sizeof(double)); 

  setupProposals(initParams); 
  
  initParamDump(); 

  for(int j = 0; j < traln->getNumberOfPartitions(); ++j)
    traln->initRevMat(j);

  evaluateFullNoBackup(this);   
  

  tree *tr = traln->getTr();
  Randomness *rand = getChainRand(); 
  stringstream ss; 
  ss << *rand ;	  	  
  printInfo("init lnl=%f\tTL=%f\tseeds=%s\n", traln->getTr()->likelihood, branchLengthToReal(tr, getTreeLength(traln, tr->nodep[1]->back)), ss.str().c_str()); 

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


/**
   @brief Normalizes the weights of the proposals in this category
 */
void Chain::normalizePropSubCats()
{
  double* catWeights = (double*)exa_calloc(NUM_PROP_CATS + 1,sizeof(double)); 

  for(int i = 0; i < this->numProposals; ++i)
    {
      proposalFunction *pf = this->proposals[i]; 
      assert(pf->category); 
      catWeights[pf->category] +=  pf->currentWeight;       
    }

  for(int i= 0; i < this->numProposals; ++i)
    {
      proposalFunction *pf = this->proposals[i]; 
      pf->currentWeight /= catWeights[pf->category]; 
    }

  exa_free(catWeights); 
}


/**
   @brief normalizes the categories  
 */
void Chain::normalizeCategories()
{
  double sum = 0.; 
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    sum += this->categoryWeights[i]; 
  for(int i = 0; i < NUM_PROP_CATS; ++i)
    this->categoryWeights[i] /= sum;   
}





/**
   @brief   Initializes the proposals based on weights given in the config file. 

   Also normalizes all weights to 1.   

 */
void Chain::setupProposals( initParamStruct *initParams)
{
  int ctr = 0; 

  proposalFunction **pfs = (proposalFunction**)exa_calloc(NUM_PROPOSALS, sizeof(proposalFunction*));  
  for(int i = 0; i < NUM_PROPOSALS; ++i)
    {
      proposalFunction *pf = NULL; 
      initProposalFunction((proposal_type)i, initParams, &pf); 
      if(pf != (proposalFunction*)NULL)
	pfs[ctr++] = pf; 
    }  
  this->numProposals = ctr; 
  this->proposals = pfs; 

  for(int i = 0; i < this->numProposals; ++i)
    {
      proposalFunction *pf = this->proposals[i]; 
      this->categoryWeights[pf->category-1] += pf->currentWeight; 
    }
  normalizeCategories();  
  normalizePropSubCats(); 
  
  // meh =/ 
  if(couplingId == 0 && runid == 0)
    printAllProposalWeights();
}



void Chain::printAllProposalWeights()
{
  if(not isOutputProcess())
    return; 
  
  printf("cat weights: TOPO=%f\tBL=%f\tFREQ=%f\tSUBST=%f\tHET=%f\n", 
	 this->categoryWeights[0],
	 this->categoryWeights[1],
	 this->categoryWeights[2],
	 this->categoryWeights[3],
	 this->categoryWeights[4] ); 

  printf("rel. prop weihgts: "); 
  for(int i = 0; i < this->numProposals; ++i)
    {
      proposalFunction *pf = this->proposals[i]; 
      printf("\t%s=%f", pf->name, pf->currentWeight);       
    }
  printf("\n"); 
}



Chain& Chain::operator=(Chain& rhs)
{
  TreeAln &thisTraln = *traln; 
  TreeAln &rhsTraln = *(rhs.traln); 

  thisTraln = rhsTraln; 

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





void Chain::debug_printAccRejc(proposalFunction *pf, bool accepted) 
{
#ifdef DEBUG_SHOW_EACH_PROPOSAL
  if(isOutputProcess())
    {
      if(accepted)
	printInfo( "ACC\t");   
      else 
	printInfo("rej\t");   	  
      printf("%s %g\n" ,pf->name, this->traln->getTr()->likelihood); 
    }
#endif
}




void Chain::drawProposalFunction(proposalFunction **result )
{    
  *result = NULL; 
  category_t
    cat = category_t(chainRand->drawSampleProportionally(categoryWeights, NUM_PROP_CATS) + 1); /* it is 1-based */
  
  double sum = 0; 
  for(int i = 0; i < numProposals; ++i)
    {
      proposalFunction *pf = proposals[i]; 
      if(pf->category == cat)
	sum += pf->currentWeight; 
    }
  assert(fabs(sum - 1.) < 0.000001); 

  double r = chainRand->drawRandDouble01();
  for(int i = 0; i < numProposals; ++i)
    {
      proposalFunction *pf = proposals[i]; 
      if(pf->category == cat)
	{
	  if(  r < pf->currentWeight)
	    {
	      *result =  pf; 
	      return; 
	    }
	  else 
	    {
	      r -= pf->currentWeight; 
	    }
	}
    }

  assert(result != NULL); 
}




void Chain::step()
{
  this->restorer->resetRestorer(*traln);

  // inform the rng that we produce random numbers for generation x  
  chainRand->rebase(currentGeneration);

  tree *tr = this->traln->getTr();   

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = getChainHeat();

  proposalFunction *pf = NULL;   
  drawProposalFunction(&pf);

  /* reset proposal ratio  */
  this->hastings = 1; 

  double oldPrior = this->priorProb; 		/* TODO  */

  /* chooses move, sets proposal ratio, correctly modifies the prior */
  pf->apply_func(this, pf);  
  double priorRatio  = this->priorProb - oldPrior; 
  /* enable once we actually have priors  */
  assert(priorRatio == 0); 

  /* chooses the cheapest way to evaluate the likelihood  */
  pf->eval_lnl(this, pf); 

  double testr = chainRand->drawRandDouble01();
  double acceptance = 
    exp((priorRatio  + tr->likelihood - prevLnl) * myHeat) 
    * this->hastings ; 

  bool wasAccepted  = testr < acceptance; 
  debug_printAccRejc( pf, wasAccepted); 

  if(wasAccepted)
    {
      pf->sCtr.accept();
      expensiveVerify(this);
    }
  else
    {
      pf->reset_func(this, pf); 
      pf->sCtr.reject();
      this->restorer->restoreArrays(*traln); // restores the previous tree state 
    }

  debug_checkTreeConsistency(this->traln->getTr());

  if( this->couplingId == 0	/* must be the cold chain  */
       && (this->currentGeneration % gAInfo.samplingFrequency) == gAInfo.samplingFrequency - 1  ) 
    {

      if( isOutputProcess() ) 
	printSample(this);       

      // TODO keep? 
      // gAInfo.bipHash->addBipartitionsToHash(*(this->traln),  this->id / gAInfo.numberCoupledChains); 
    }


  /* autotuning for proposal parameters. With increased parallelism
     this will become more complicated.  */
  if( gAInfo.tuneFreq > 0 && this->currentGeneration % gAInfo.tuneFreq == gAInfo.tuneFreq - 1 )
    {
      for(int i = 0; i < this->numProposals; ++i)
	{
	  proposalFunction *pf = this->proposals[i]; 
	  
	  if(pf->autotune)	/* only, if we set this   */
	    {
#ifdef TUNE_ONLY_IF_ENOUGH
	      if(pf->sCtr.lAcc + pf->sCtr.lRej < gAInfo.tuneFreq )
		continue; 
#endif
	      pf->autotune(this, pf);
	    }
	}
    }

  this->currentGeneration++; 
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


