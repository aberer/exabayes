#include "Chain.hpp"
#include "Topology.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp"
#include "randomness.h"
// #include "chain.h"
#include "output.h"
#include "globals.h"
#include "BipartitionHash.hpp"
#include "adapters.h"
#include "proposals.h"
#include "topology-utils.h"



// TODO 
#define TOPO_RESTORE 0 
#define TOPO_SAVE 1 





Chain::Chain(randKey_t seed, int id, int runid, TreeAln* _traln, initParamStruct *initParams)  
  : traln(_traln)
  , couplingId(id)
  , currentGeneration(0)
  , hastings(1)
{
  rKey.v[0] = seed.v[0]; 
  rKey.v[1] = seed.v[1]; 

  // make the file names 
  
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
  saveTreeStateToChain(); 
  PRINT("Initial LnL for chain %d is  %f\ttree-length=%.3f\tseed=%u,%u\n", this->id, this->traln->getTr()->likelihood, 
	branchLengthToReal(this->traln->getTr(), getTreeLength(this->traln, this->traln->getTr()->nodep[1]->back)),
	this->rKey.v[0], this->rKey.v[1]
	); 

  // TODO  
  assert(0); 
}



/**
   @brief returns the inverse temperature for this chain
 */ 
double Chain::getChainHeat()
{
  double tmp = 1. + deltaT * couplingId; 
  double inverseHeat = 1 / tmp; 
  assert(inverseHeat < 1.); 
  return inverseHeat; 
}




// LEGACY
static void traverseAndTreatBL(node *p, TreeAln *traln, double *blBuf, int* cnt, boolean restore)
{
  tree *tr = traln->getTr();
  nodeptr q; 
  assert(traln->getNumBranches() == 1); 

  if(restore == TOPO_RESTORE )
    traln->setBranchLengthSave(blBuf[*cnt], 0,p); 
  else if(restore == TOPO_SAVE)    
    blBuf[*cnt] =  traln->getBranchLength(p,0); 
  else 
    assert(0); 
  *cnt += 1 ; 
  
  if( NOT isTip(p->number, tr->mxtips))
    {
      q = p->next; 
      while(q != p )
	{
	  traverseAndTreatBL(q->back, traln, blBuf, cnt, restore); 
	  q = q->next; 
	}
    }
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
  
  /* only print that once  */
  if(this->id == 0)
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




/**
   @brief Save all relevan information from the tree into the chain
   Chain. 
 */ 
void Chain::saveTreeStateToChain()
{
  tree *tr  = this->traln->getTr();     
  dump.topology->saveTopology(*(traln));

  /* save branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, this->traln, this->dump.branchLengths, &cnt, TOPO_SAVE); 
  assert(cnt == 2 * tr->mxtips - 3 ); 

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
  tree *tr = traln->getTr(); 
  
  /* TODO enable multi-branch    */
  assert(traln->getNumBranches() == 1); 

  dump.topology->restoreTopology(*(traln));

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, traln, dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 

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

  wasAccepted = TRUE; 

}






void Chain::drawProposalFunction(proposalFunction **result )
{    
  *result = NULL; 
  category_t
    cat = category_t(drawSampleProportionally(this,categoryWeights, NUM_PROP_CATS) + 1); /* it is 1-based */

  /* printInfo(chain, "drawing proposal; category is %d\n"), cat;  */
  
  double sum = 0; 
  for(int i = 0; i < numProposals; ++i)
    {
      proposalFunction *pf = proposals[i]; 
      if(pf->category == cat)
	sum += pf->currentWeight; 
    }
  assert(fabs(sum - 1.) < 0.000001); 

  double r = drawRandDouble01(this);
  /* printf("numProp=%d\n", chain->numProposals);  */
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
  this->restorer->resetRestorer();

  tree *tr = this->traln->getTr();   

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = this->getChainHeat();

  // double myHeat = getChainHeat(chain ) ; 

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

  double testr = drawRandDouble01(this);
  double acceptance = 
    exp((priorRatio  + tr->likelihood - prevLnl) * myHeat) 
    * this->hastings ; 

  this->wasAccepted  = testr < acceptance; 
  debug_printAccRejc(this, pf, this->wasAccepted); 

  if(this->wasAccepted)
    {
      pf->sCtr.accept();
      expensiveVerify(this);
    }
  else
    {
      pf->reset_func(this, pf); 
      pf->sCtr.reject();
      this->restorer->restore(); // restores the previous tree state 

    }

  debug_checkTreeConsistency(this->traln->getTr());

  if( this->couplingId == 0	/* must be the cold chain  */
       && (this->currentGeneration % gAInfo.samplingFrequency) == gAInfo.samplingFrequency - 1  ) 
    {

      if( isOutputProcess() ) 
	printSample(this);       

      // TODO keep? 
      gAInfo.bipHash->addBipartitionsToHash(*(this->traln),  this->id / gAInfo.numberCoupledChains); 
    }


  /* the output for the console  */
  if(isOutputProcess() 
     && this->couplingId == 0     
     && gAInfo.printFreq > 0 
     && this->currentGeneration % gAInfo.printFreq == gAInfo.printFreq -1   )
    {
      chainInfo(this); 
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





void Chain::initParamDump( )
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

  dump.branchLengths = (double*)exa_calloc(2 * tr->mxtips, sizeof(double));
  for(int i = 0; i < 2 * tr->mxtips; ++i)
    dump.branchLengths[i] = TreeAln::initBL;   
}
