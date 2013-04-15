
#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "main-common.h"
#include "chain.h"
#include "proposals.h"
#include "output.h"
#include "convergence.h" 
#include "treeRead.h"
#include "topology-utils.h"
#include "randomness.h"
#include "eval.h"
#include "adapters.h"
#include "randomTree.h" 

#include "TreeAln.hpp"
#include "LnlRestorer.hpp"
#include "Topology.hpp"


void initDefaultValues(state *chain, tree *tr)
{
  /* TODO we have to redo the prior framework */
  chain->priorProb = 1; 
  chain->hastings = 1; 
  chain->currentGeneration = 0; 

  chain->penaltyFactor = 0.0;
}


void makeChainFileNames(state *chain, int num)
{
  char tName[1024],
    pName[1024] ; 
  
  sprintf(tName, "%s%s_topologies.%s.%d", workdir, PROGRAM_NAME, run_id, num); 
  sprintf(pName, "%s%s_parameters.%s.%d", workdir, PROGRAM_NAME, run_id, num); 
  
  /* todo binary state file?  */
  chain->topologyFile = fopen(tName, "w"); 
  chain->outputParamFile = fopen(pName, "w");   
}



/**
   @brief returns the inverse temperature for this chain
 */ 
double getChainHeat(state *chain )
{
  int runId = chain->id / gAInfo.numberCoupledChains;  
  const double deltaT = gAInfo.temperature[runId]; 

  if(chain->couplingId == 0 )
    return 1; 
  
  double tmp  = 1. + deltaT * chain->couplingId; 
  double myHeat = 1. / (double)tmp; 

  assert(myHeat < 1.);
  return myHeat; 
}


/**
   @brief corrects branch lengths when those are read from a file
*/
void traverseInitCorrect(nodeptr p, int *count, TreeAln *traln )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;

  for( i = 0; i < traln->getNumBranches(); i++)
    p->z[i] =  exp( - p->z[i] / tr->fracchange); 
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitCorrect(q->back, count, traln);
	  q = q->next;
	} 
    }
}


static void initParamDump(TreeAln *traln, paramDump *dmp)
{  
  tree *tr = traln->getTr();
  
  dmp->topology = new Topology(tr->mxtips); 
  dmp->infoPerPart = (perPartitionInfo*)exa_calloc(traln->getNumberOfPartitions(), sizeof(perPartitionInfo));
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *p =  dmp->infoPerPart + i ; 
      p->alpha = 0.5; 
      for(int j = 0; j < 5; ++j)
	p->substRates[j] = 0.5; 
      p->substRates[5] = 1; 
      for(int j = 0; j < 4; ++j)
	p->frequencies[j] = 0.25;       
    }

  dmp->branchLengths = (double*)exa_calloc(2 * tr->mxtips, sizeof(double));
  for(int i = 0; i < 2 * tr->mxtips; ++i)
    dmp->branchLengths[i] = INIT_BRANCH_LENGTHS;   
}




static void copyParamDump(TreeAln *traln,paramDump *dest, const paramDump *src)
{
  /* no topo or branch lengths! this should be done by copytopo */
  
  for(int i = 0; i < traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *pDest =  dest->infoPerPart + i ,
	*pSrc = src->infoPerPart + i; 
      
      pDest->alpha = pSrc->alpha; 
      
      for(int j = 0; j < 6; ++j)
	pDest->substRates[j] = pSrc->substRates[j]; 
      for(int j = 0; j < 4; ++j)
	pDest->frequencies[j] = pSrc->frequencies[j]; 
    }

}




void copyState(state *dest, const state *src )
{
  TreeAln &thisTraln =  *(dest->traln); 
  TreeAln &rhsTraln = *(src->traln); 

  thisTraln = rhsTraln; 
  *(dest->traln ) = *(src->traln); 
  
  dest->dump.topology->saveTopology(thisTraln);
}




static void setupGlobals(initParamStruct *initParams, tree *tr)
{
  if (initParams->numGen > 0)
    gAInfo.numGen = initParams->numGen; 

  gAInfo.samplingFrequency = initParams->samplingFrequency; 
  gAInfo.diagFreq = initParams->diagFreq; 
  gAInfo.numberOfRuns =   initParams->numIndiChains; 
  gAInfo.numberCoupledChains = initParams->numCoupledChains; 



  /* initialize a matrix of swaps (wasting some space here) */
  gAInfo.swapInfo = (successCtr**)exa_calloc(gAInfo.numberOfRuns, sizeof(successCtr*)); 
  int n = gAInfo.numberCoupledChains; 
  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    gAInfo.swapInfo[i] = (successCtr*)exa_calloc( n * n , sizeof(successCtr)); 

  gAInfo.temperature = (double*)exa_calloc(gAInfo.numberOfRuns, sizeof(double)); 
  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    gAInfo.temperature[i] = HEAT_FACTOR; 

  unsigned int bvLength = 0; 
  tr->bitVectors = initBitVector(tr->mxtips, &bvLength); 
  hashtable *ht = initHashTable(tr->mxtips * tr->mxtips * 10);
  gAInfo.bvHash = ht; 
}


#include <vector>
using namespace std; 


/**
   @brief An overloaded initialization function for all the chains. 

   This function should initialize all chain data structures and parse
   the config file.

   This function also decides which aln,tr structures are assigned to
   which chains.
 */ 
void initializeIndependentChains(tree *tr, analdef *adef, state **resultIndiChains)
{
  FILE *treeFH = NULL; 
  if( gAInfo.numberOfStartingTrees > 0 )
    treeFH = myfopen(tree_file, "r"); 

  initParamStruct *initParams = (initParamStruct*)exa_calloc(1,sizeof(initParamStruct));   
  parseConfigWithNcl(configFileName, &initParams);  
  setupGlobals(initParams, tr); 
  int totalNumChains = gAInfo.numberOfRuns * gAInfo.numberCoupledChains;   

  PRINT("number of independent runs=%d, number of coupled chains per run=%d => total of %d chains \n", gAInfo.numberOfRuns, gAInfo.numberCoupledChains, totalNumChains ); 
  *resultIndiChains = (state*)exa_calloc( totalNumChains , sizeof(state));     

#ifdef DEBUG_LNL_VERIFY
  gAInfo.debugTree = new TreeAln(tr ); 
#endif

  for(int i = 0; i < totalNumChains; ++i)
    {       
      state *theChain = *resultIndiChains + i;       
      if( i < gAInfo.numberCoupledChains)
	theChain->traln = new TreeAln(tr); 
      else 
	{
	  state *masterChain = *resultIndiChains + (i %  gAInfo.numberCoupledChains) ; 
	  theChain->traln = new TreeAln(*(masterChain->traln)); // initialize from master chain 
	}
      
      TreeAln *traln = 	theChain->traln; 
      theChain->id = i; 
      theChain->couplingId = i % gAInfo.numberCoupledChains ; 

      theChain->categoryWeights = (double*)exa_calloc(NUM_PROP_CATS, sizeof(double)); 

      tree *myTree = theChain->traln->getTr();       
      initDefaultValues(theChain, myTree);

      if(theChain->id / gAInfo.numberCoupledChains == 0 ) 
	theChain->restorer = new LnlRestorer(theChain);
      else 
	{	  
	  state *leechChain =   (*resultIndiChains ) + theChain->couplingId; 
	  theChain->restorer = leechChain->restorer; 
	}

      setupProposals(theChain, initParams); 

      /* init the param dump  */      
      initParamDump(traln, &(theChain->dump)); 

      for(int j = 0; j < traln->getNumberOfPartitions(); ++j ) 
	theChain->traln->initRevMat(j);

      initLocalRng(theChain); 
      
      if( i % gAInfo.numberCoupledChains == 0)
	{
	  /* initialize with new tree  */
	  if( i  / gAInfo.numberCoupledChains < gAInfo.numberOfStartingTrees )
	    {
	      boolean hasBranchLength = readTreeWithOrWithoutBL(theChain->traln->getTr(), treeFH); 

	      if(hasBranchLength)
		{
		  int count = 0; 
		  traverseInitCorrect(tr->start->back, &count, theChain->traln ) ;
		  assert(count == 2 * theChain->traln->getTr()->mxtips -3); 
		}
	      else 
		{
		  /* set some standard branch lengths */
		  /* TODO based on prior?   */
		  int count = 0; 
		  traverseInitFixedBL( theChain->traln->getTr()->start->back, &count, theChain->traln, INIT_BRANCH_LENGTHS );
		  assert(count  == 2 * theChain->traln->getTr()->mxtips - 3);	  	      
		}

	      if(processID == 0)
		printf("initializing chain %d with provided starting tree\n", i); 
	    }
	  else 
	    {
	      exa_makeRandomTree(theChain->traln->getTr());
	      if(processID == 0)
		printf("initializing chain %d with random tree\n", i); 

	      /* TODO maybe prior for initial branch lengths   */
	      int count = 0; 
	      traverseInitFixedBL( theChain->traln->getTr()->start->back, &count, theChain->traln, INIT_BRANCH_LENGTHS );
	      assert(count  == 2 * tr->mxtips - 3);	  
	    }
	}
      else 
	{
	  state *coldChain = *resultIndiChains + (theChain->id / gAInfo.numberCoupledChains) * gAInfo.numberCoupledChains; 
	  copyState(theChain , coldChain);
	  applyChainStateToTree(theChain); 
	}

      evaluateFullNoBackup(theChain);

      /* now save the tree to a chain chains */
      saveTreeStateToChain(theChain); 

      PRINT("init lnl for chain %d is  %f\ttree-length=%.3f\n", theChain->id, theChain->traln->getTr()->likelihood, 
	    branchLengthToReal(theChain->traln->getTr(), getTreeLength(theChain->traln->getTr(), theChain->traln->getTr()->nodep[1]->back))); 

      if(processID == 0 )
	{	  
	  if( i % gAInfo.numberCoupledChains == 0)
	    {
	      makeChainFileNames(theChain, i / gAInfo.numberCoupledChains); 
	      initializeOutputFiles(theChain); 
	    }
	  else 
	    {	      
	      state *previousChain = *resultIndiChains + (i  - 1 ) ; 
	      theChain->topologyFile = previousChain->topologyFile; 
	      theChain->outputParamFile = previousChain->topologyFile; 
	    }
	}
    }

  if(gAInfo.numberOfStartingTrees > 0)
    fclose(treeFH); 

  // TODO kill initial tree 
}



void traverseInitFixedBL(nodeptr p, int *count, TreeAln *traln,  double z )
{
  tree *tr = traln->getTr();
  nodeptr q;
  int i;
  
  for( i = 0; i < traln->getNumBranches(); i++)
    p->z[i] = p->back->z[i] = z;  
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitFixedBL(q->back, count, traln, z);
	  q = q->next;
	} 
    }
}



void traverseAndPrint(nodeptr p, int *count, tree *tr)
{
  nodeptr q; 
  printf("%d:%f,", p->number, p->z[0]); 

  if( ! isTip(p->number, tr->mxtips))
    {
      q = p->next; 
      while(p != q)
	{
	  traverseAndPrint(q->back,count, tr); 
	  q = q->next; 
	}
    }  
}



void traverseAndTreatBL(node *p, TreeAln *traln, double *blBuf, int* cnt, boolean restore)
{
  tree *tr = traln->getTr();
  nodeptr q; 
  assert(traln->getNumBranches() == 1); 

  if(restore == TOPO_RESTORE )
    { 
      /* TODO kassian: is this correct? */
      p->z[0] = blBuf[*cnt];             
      p->back->z[0] = blBuf[*cnt]; 
      
    }
  else if(restore == TOPO_SAVE)    
    {
      blBuf[*cnt] =  p->z[0]; 
    }
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
   @brief Applies the state of the chain to its tree. 

   Notice: you must not simply change the tree pointer. More
   modifications are necessary to do so.

   @param boolean checkLnl -- should we check, if the lnl is the same as
   before? If we applied it the first time, there is no before.
 */ 
void applyChainStateToTree(state *chain)
{
  tree *tr = chain->traln->getTr(); 
  
  /* TODO enable multi-branch    */
  assert(chain->traln->getNumBranches() == 1); 

  chain->dump.topology->restoreTopology(*(chain->traln));

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, chain->traln, chain->dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 

  /* restore model parameters */
  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i ; 

      pInfo *partition = chain->traln->getPartition( i);
      partition->alpha = info->alpha;

      memcpy(partition->substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(partition->frequencies, info->frequencies, 4 * sizeof(double));
      
      chain->traln->initRevMat(i);
      makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);
    } 

  // evaluateGenericWrapper(chain, tr->start, TRUE ); 
  evaluateFullNoBackup(chain); 
  
  if(processID == 0)
    {
#ifdef DEBUG_BL
      printf("branch lengths: ");
      int count = 0; 
      traverseAndPrint(tr->start->back, &count, tr); 
      printf("\n"); 
#endif
    }


  /* TODO for now  */
  chain->wasAccepted = TRUE; 
  // chain->prevProposal = NULL; 

#if 0 				// TODO 
  updateSavedState(chain); 
#endif
}



/**
   @brief Save all relevan information from the tree into the chain
   state. 
 */ 
void saveTreeStateToChain(state *chain)
{
  tree *tr  = chain->traln->getTr();     
  chain->dump.topology->saveTopology(*(chain->traln));

  /* save branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, chain->traln, chain->dump.branchLengths, &cnt, TOPO_SAVE); 
  assert(cnt == 2 * tr->mxtips - 3 ); 

  /* save model parameters */
  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i; 
      pInfo *partition = chain->traln->getPartition(i); 
      
      info->alpha =  partition->alpha; 
      
      memcpy(info->substRates, partition->substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, partition->frequencies, 4 * sizeof(double)); 
    }  


  if(processID == 0)
    {

#ifdef DEBUG_BL
      int count = 0;       
      printf("saving state: ");
      traverseAndPrint(tr->start->back,  &count,tr);
      printf("\n");
#endif
    }
}





/**
   @brief draws a proposal function.

   Notice: this could be extended later, if we decide to make this
   dependent on the previous state.
   
   Furthermore, we must be sure now that category weights and relative
   proposal weights sum up to 1 each. 
   
 */ 
void drawProposalFunction(state *chain, proposalFunction **result )
{  
  
  *result = NULL; 
  category_t
    cat = category_t(drawSampleProportionally(chain,chain->categoryWeights, NUM_PROP_CATS) + 1); /* it is 1-based */

  /* printInfo(chain, "drawing proposal; category is %d\n"), cat;  */
  
  double sum = 0; 
  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
      if(pf->category == cat)
	sum += pf->currentWeight; 
    }
  assert(fabs(sum - 1.) < 0.000001); 

  double r = drawRandDouble01(chain);
  /* printf("numProp=%d\n", chain->numProposals);  */
  for(int i = 0; i < chain->numProposals; ++i)
    {
      proposalFunction *pf = chain->proposals[i]; 
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


/**
   @brief Execute one generation of a given chain.  
 */
void step(state *chain)
{
  chain->restorer->resetRestorer();

  tree *tr = chain->traln->getTr();   

  assert(tr->fracchange > 0); 

  double prevLnl = tr->likelihood;     

  double myHeat = getChainHeat(chain ) ; 

  proposalFunction *pf = NULL;   
  drawProposalFunction(chain, &pf);

  /* reset proposal ratio  */
  chain->hastings = 1; 

  double oldPrior = chain->priorProb; 		/* TODO  */

  /* chooses move, sets proposal ratio, correctly modifies the prior */
  pf->apply_func(chain, pf);  
  double priorRatio  = chain->priorProb - oldPrior; 
  /* enable once we actually have priors  */
  assert(priorRatio == 0); 

  /* chooses the cheapest way to evaluate the likelihood  */
  pf->eval_lnl(chain, pf); 

  double testr = drawRandDouble01(chain);
  double acceptance = 
    exp((priorRatio  + tr->likelihood - prevLnl) * myHeat) 
    * chain->hastings ; 

  chain->wasAccepted  = testr < acceptance; 
  debug_printAccRejc(chain, pf, chain->wasAccepted); 

  if(chain->wasAccepted)
    {
      cntAccept(&(pf->sCtr)); 
      expensiveVerify(chain);
    }
  else
    {
      pf->reset_func(chain, pf); 
      cntReject(&(pf->sCtr));       
      chain->restorer->restore(); // restores the previous tree state 

    }

  debug_checkTreeConsistency(chain->traln->getTr());

  if(  processID == 0 
       && chain->couplingId == 0	/* must be the cold chain  */
       && (chain->currentGeneration % gAInfo.samplingFrequency) == gAInfo.samplingFrequency - 1  ) 
    {
      printSample(chain);       

      if(chain->currentGeneration >  BURNIN)
	addBipartitionsToHash(tr, chain ); 
    }

  
  /* the output for the console  */
  if(processID == 0 
     && chain->couplingId == 0     
     && chain->currentGeneration % PRINT_FREQUENCY == PRINT_FREQUENCY -1   )
    {
      chainInfo(chain); 
    }


#ifdef TUNE_PARAMETERS
  /* autotuning for proposal parameters. With increased parallelism
     this will become more complicated.  */
  if( chain->currentGeneration % TUNE_FREQUENCY == TUNE_FREQUENCY - 1 )
    {
      for(int i = 0; i < chain->numProposals; ++i)
	{
	  proposalFunction *pf = chain->proposals[i]; 
	  
	  if(pf->autotune)	/* only, if we set this   */
	    {
#ifdef TUNE_ONLY_IF_ENOUGH
	      if(pf->sCtr.lAcc + pf->sCtr.lRej < TUNE_FREQUENCY )
		continue; 
#endif
	      pf->autotune(chain, pf);
	    }
	}

    }
#endif

  chain->currentGeneration++; 
}
