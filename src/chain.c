


#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "main-common.h"
#include "chain.h"
#include "proposals.h"
#include "output.h"
#include "convergence.h" 
#include "treeRead.h"
#include "exa-topology.h"
#include "randomness.h"
#include "eval.h"
#include "adapters.h"
#include "randomTree.h" 




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
void traverseInitCorrect(nodeptr p, int *count, tree *tr )
{
  nodeptr q;
  int i;
  
  assert(NOT hasPergeneBL(tr)); 
  
  for( i = 0; i < getNumBranches(tr); i++)
    p->z[i] =  exp( - p->z[i] / tr->fracchange); 
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitCorrect(q->back, count, tr);
	  q = q->next;
	} 
    }
}


static void initParamDump(tree *tr, paramDump *dmp)
{  
  dmp->topo = setupTopol(tr->mxtips); 
  dmp->infoPerPart = exa_calloc(getNumberOfPartitions(tr), sizeof(perPartitionInfo));
  for(int i = 0; i < getNumberOfPartitions(tr); ++i)
    {
      perPartitionInfo *p =  dmp->infoPerPart + i ; 
      p->alpha = 0.5; 
      for(int j = 0; j < 5; ++j)
	p->substRates[j] = 0.5; 
      p->substRates[5] = 1; 
      for(int j = 0; j < 4; ++j)
	p->frequencies[j] = 0.25;       
    }

  dmp->branchLengths = exa_calloc(2 * tr->mxtips, sizeof(double));
  for(int i = 0; i < 2 * tr->mxtips; ++i)
    dmp->branchLengths[i] = INIT_BRANCH_LENGTHS;   
}




static void copyParamDump(tree *tr,paramDump *dest, const paramDump *src)
{
  /* no topo or branch lengths! this should be done by copytopo */
  
  for(int i = 0; i < getNumberOfPartitions(tr); ++i)
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
  copyTopology(dest->tr, src->tr); 
  copyParamDump(dest->tr, &(dest->dump), &(src->dump)); 
  saveTree(dest->tr, dest->dump.topo); 	  
  /* applyChainStateToTree(dest);  */
}



void preinitTree(tree *tr)
{   
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA; 
  tr->multiStateModel  = GTR_MULTI_STATE;
#if (HAVE_PLL == 0 ) 
    tr->useGappedImplementation = FALSE;
    tr->saveBestTrees          = 0;
#endif
  tr->saveMemory = FALSE;
  tr->manyPartitions = FALSE;
  tr->categories             = 25;
  tr->grouped = FALSE;
  tr->constrained = FALSE;
  tr->gapyness               = 0.0; 
  tr->useMedian = FALSE;
}



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

  initParamStruct *initParams = exa_calloc(1,sizeof(initParamStruct)); 
  
  parseConfigWithNcl(configFileName, &initParams);
  
  if (initParams->numGen > 0)
    gAInfo.numGen = initParams->numGen; 

  gAInfo.samplingFrequency = initParams->samplingFrequency; 
  gAInfo.diagFreq = initParams->diagFreq; 
  gAInfo.numberOfRuns =   initParams->numIndiChains; 
  gAInfo.numberCoupledChains = initParams->numCoupledChains; 
  int totalNumChains = gAInfo.numberOfRuns * gAInfo.numberCoupledChains; 

  /* initialize a matrix of swaps (wasting some space here) */
  gAInfo.swapInfo = (successCtr**)exa_calloc(gAInfo.numberOfRuns, sizeof(successCtr*)); 
  int n = gAInfo.numberCoupledChains; 
  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    gAInfo.swapInfo[i] = exa_calloc( n * n , sizeof(successCtr)); 

  gAInfo.temperature = exa_calloc(gAInfo.numberOfRuns, sizeof(double)); 
  for(int i = 0; i < gAInfo.numberOfRuns; ++i)
    gAInfo.temperature[i] = HEAT_FACTOR; 


#ifdef MC3_SPACE_FOR_TIME
  int treesNeeded = gAInfo.numberCoupledChains  -1 ; 
#if HAVE_PLL == 1 
  gAInfo.partitions = (partitionList*)exa_realloc(gAInfo.partitions, (treesNeeded + 1)  * sizeof(partitionList)); 
#endif
  tree *moreTrees = exa_calloc(treesNeeded, sizeof(tree)); 
  for(int i = 0; i < treesNeeded; ++i)
    {
      preinitTree(moreTrees+i); 
#if HAVE_PLL == 1 
      initializeTree(moreTrees + i, gAInfo.partitions + i+1 ,adef); 
#else 
      initializeTree(moreTrees + i ,adef); 
#endif
      tree *oT = moreTrees + i; 
      oT->mxtips = tr->mxtips; 
    }
  
  for(int i = 0; i < treesNeeded; ++i)
    {
      tree *trH = moreTrees + i ; 
      for(int j = 1; j <= trH->mxtips; ++j)
	trH->nodep[j]->hash = tr->nodep[j]->hash; 
    }
#endif
  
  if(processID == 0)
    printf("number of independent runs=%d, number of coupled chains per run=%d => total of %d chains \n", gAInfo.numberOfRuns, gAInfo.numberCoupledChains, totalNumChains ); 
  *resultIndiChains = exa_calloc( totalNumChains , sizeof(state));     
  
  unsigned int bvLength = 0; 
  tr->bitVectors = initBitVector(tr->mxtips, &bvLength); 

  hashtable *ht = initHashTable(tr->mxtips * tr->mxtips * 10);
  gAInfo.bvHash = ht; 
  int numPart = getNumberOfPartitions(tr); 
  
  for(int i = 0; i < totalNumChains; ++i)
    { 
      state *theChain = *resultIndiChains + i;     
      
      theChain->id = i; 
      theChain->couplingId = i % gAInfo.numberCoupledChains ; 

      theChain->categoryWeights = exa_calloc(NUM_PROP_CATS, sizeof(double)); 
      
#ifdef MC3_SPACE_FOR_TIME
      /* important NOTICE : in this scheme, we assume, that there is
	 one tree (partitionList) for each coupled chain (!= total
	 number of chains) */
      assert(theChain->couplingId < gAInfo.numberCoupledChains); 
      theChain->tr = (theChain->couplingId ==  0) ? tr : moreTrees +  theChain->couplingId - 1 ; 
#else 
      /* here we are only working with one single tree */	      
      theChain->tr = tr;

#endif 

#if HAVE_PLL == 1 
#ifdef MC3_SPACE_FOR_TIME
      theChain->partitions = gAInfo.partitions + theChain->couplingId ; 
#else 
      theChain->partitions = gAInfo.partitions; 
#endif
#endif

      tree *myTree = theChain->tr; 
      myTree->bitVectors = tr->bitVectors; 
      
      initDefaultValues(theChain, myTree);

      if(theChain->id / gAInfo.numberCoupledChains == 0 ) 
	{
	  theChain->lnl.orientation = exa_calloc(tr->mxtips,sizeof(int)); 
	  theChain->lnl.vectorsPerPartition = exa_calloc(numPart, sizeof(double**)); 
	  theChain->lnl.partitionScaler = (nat**)exa_calloc(numPart, sizeof(nat*)); 	  

	  for(int j = 0; j < numPart; ++j)
	    {
	      pInfo *partition = getPartition(theChain, j); 
	      theChain->lnl.vectorsPerPartition[j] = exa_calloc(tr->mxtips, sizeof(double*)); 
#if HAVE_PLL == 1 		/* TODO that's hacky */
	      int length = partition->upper - partition->lower; 
#else 
	      int length = partition->width; 
#endif
	      /* printf(" [%d] length is %d (width=%d)\n", processID, length, (int)partition->width) ; */
	      assert(length > 0); 
	      for(int k = 0; k < tr->mxtips; ++k) 
		theChain->lnl.vectorsPerPartition[j][k] = exa_calloc( length *  LENGTH_LNL_ARRAY , sizeof(double)); 

	      theChain->lnl.partitionScaler[j] = (nat*)exa_calloc(2 * tr->mxtips,sizeof(nat)); 	      
	    }	  
	}
      else 
	{
	  state *leechChain =   (*resultIndiChains ) + theChain->couplingId; 
	  theChain->lnl.orientation = leechChain->lnl.orientation; 
	  theChain->lnl.vectorsPerPartition = leechChain->lnl.vectorsPerPartition; 
	  theChain->lnl.partitionScaler = leechChain->lnl.partitionScaler; 
	}
      
      setupProposals(theChain, initParams); 

      /* init the param dump  */
      initParamDump(myTree, &(theChain->dump)); 

      for(int i = 0; i < getNumberOfPartitions(myTree); ++i ) 
	exa_initReversibleGTR(theChain,i);

      initLocalRng(theChain); 
      
      if( i % gAInfo.numberCoupledChains == 0)
	{
	  /* initialize with new tree  */
	  if( i  / gAInfo.numberCoupledChains < gAInfo.numberOfStartingTrees )
	    {
	      boolean hasBranchLength = readTreeWithOrWithoutBL(theChain->tr, treeFH); 

	      if(hasBranchLength)
		{
		  int count = 0; 
		  traverseInitCorrect(tr->start->back, &count, theChain->tr ) ;
		  assert(count == 2 * theChain->tr->mxtips -3); 
		}
	      else 
		{
		  /* set some standard branch lengths */
		  /* TODO based on prior?   */
		  int count = 0; 
		  traverseInitFixedBL( theChain->tr->start->back, &count, theChain->tr, INIT_BRANCH_LENGTHS );
		  assert(count  == 2 * theChain->tr->mxtips - 3);	  	      
		}

	      if(processID == 0)
		printf("initializing chain %d with provided starting tree\n", i); 
	    }
	  else 
	    {
	      exa_makeRandomTree(theChain->tr);
	      if(processID == 0)
		printf("initializing chain %d with random tree\n", i); 

	      /* TODO maybe prior for initial branch lengths   */
	      int count = 0; 
	      traverseInitFixedBL( theChain->tr->start->back, &count, theChain->tr, INIT_BRANCH_LENGTHS );
	      assert(count  == 2 * tr->mxtips - 3);	  
	    }
	}
      else 
	{
#ifdef MC3_SPACE_FOR_TIME
	  state *coldChain = *resultIndiChains + (theChain->id / gAInfo.numberCoupledChains) * gAInfo.numberCoupledChains; 
	  copyState(theChain , coldChain);
	  applyChainStateToTree(theChain); 
#else 
	  /* initialize with previous state */
	  applyChainStateToTree( (*resultIndiChains) + i-1); 
#endif
	}

      evaluateGenericWrapper(theChain, theChain->tr->start, TRUE); 

      /* now save the tree to a chain chains */
      saveTreeStateToChain(theChain); 

      if(processID == 0)
	printf("init lnl for chain %d is  %f\n", theChain->id, theChain->tr->likelihood); 

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
}



void traverseInitFixedBL(nodeptr p, int *count, tree *tr,  double z )
{
  nodeptr q;
  int i;
  
  for( i = 0; i < getNumBranches(tr); i++)
    p->z[i] = p->back->z[i] = z;  
  *count += 1;
  
  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseInitFixedBL(q->back, count, tr, z);
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



void traverseAndTreatBL(node *p, tree *tr, double *blBuf, int* cnt, boolean restore)
{
  nodeptr q; 
  assert(getNumBranches(tr) == 1); 

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
	  traverseAndTreatBL(q->back, tr, blBuf, cnt, restore); 
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
  tree *tr = chain->tr; 
  
  /* TODO enable multi-branch    */
  assert(getNumBranches(tr) == 1); 

  boolean treeWasRestored = restoreTree(chain->dump.topo, tr); 
  assert(treeWasRestored);   

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 

  /* restore model parameters */
  for(int i = 0; i < getNumberOfPartitions(tr); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i ; 

      pInfo *partition = getPartition(chain, i);
      partition->alpha = info->alpha;

      memcpy(partition->substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(partition->frequencies, info->frequencies, 4 * sizeof(double));

      exa_initReversibleGTR (chain, i); 

      makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);
    } 

  evaluateGenericWrapper(chain, tr->start, TRUE ); 
  
  if(processID == 0)
    {
#ifdef DEBUG_BL
      printf("branch lengths: ");
      int count = 0; 
      traverseAndPrint(tr->start->back, &count, tr); 
      printf("\n"); 
#endif
    }


#if 0 
  if( chain->tr->likelihood != 0 && fabs (tr->likelihood - chain->lnl.likelihood ) > 1e-6 )
    {
      printInfo(chain, "WARNING: obtained a different likelihood  after restoring previous chain state (before/after): %f / %f\n", 
		chain->id, chain->currentGeneration, chain->lnl.likelihood, tr->likelihood); 
      assert( fabs(chain->lnl.likelihood - tr->likelihood ) < 1e-6 ) ; 
    }
#endif

  /* TODO for now  */
  chain->wasAccepted = TRUE; 
  chain->prevProposal = NULL; 

  saveAlignAndTreeState(chain); 
}



/**
   @brief Save all relevan information from the tree into the chain
   state. 
 */ 
void saveTreeStateToChain(state *chain)
{
  tree *tr  = chain->tr;   
  saveTree(tr, chain->dump.topo);

  /* save branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_SAVE); 
  assert(cnt == 2 * tr->mxtips - 3 ); 

  /* save model parameters */
  for(int i = 0; i < getNumberOfPartitions(tr); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i; 
      pInfo *partition = getPartition(chain,i); 
      
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
    cat = drawSampleProportionally(chain,chain->categoryWeights, NUM_PROP_CATS) + 1; /* it is 1-based */

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



void debug_checkLikelihood(state *chain)
{
#if 0 
#ifdef DEBUG_LNL_VERIFY
  tree *tr = chain->tr; 
  exa_evaluateGeneric(chain, tr->start, TRUE );   
  double diff =  fabs(chain->lnl.likelihood - tr->likelihood ); 
  if( diff > ACCEPTED_LIKELIHOOD_EPS)
    {
      printf("LNL difference: %f\n", diff  ); 
      assert(diff  < 1e-6); 
    }  
  chain->lnl.likelihood = tr->likelihood; 
#endif  
#endif
}


/**
   @brief Execute one generation of a given chain.  
 */
void step(state *chain)
{
  tree *tr = chain->tr;   

/*   if(processID == 0) */
/*     { */
/*       pInfo *partition = getPartition(chain,0 ) ;  */
/*       assert(getNumBranches(tr) == 1 );  */

/*       double indiFC = 0, */
/* 	indiFC2 = 0;  */
/* #if HAVE_PLL == 1  */
/*       indiFC = partition->fracchange;  */
/*       indiFC2 = getPartition(chain,1)->fracchange;  */
/* #else  */
/*       indiFC = tr->fracchanges[0];  */
/*       indiFC2 = tr->fracchanges[1];  */
/* #endif */

/*       printf("fracchange = %g, fracchanges[0]=%g\tfracchanges[1]=%g\tbl=%g\t", tr->fracchange, indiFC, indiFC2 , tr->nodep[1]->z[0]);  */
/*       for(int i = 0;  i < 6 ; ++i) */
/* 	printf("%g,", partition->substRates[i]);  */
/*       printf("\t");  */
/*       for(int i = 0; i < 4; ++i) */
/* 	printf("%g,", partition->frequencies[i]);  */
/*       printf("\n");  */
/*     } */
    

  assert(chain->tr->fracchange > 0); 

  double prevLnl = chain->tr->likelihood;     

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
    exp((priorRatio  + chain->tr->likelihood - prevLnl) * myHeat) 
    * chain->hastings ; 

  chain->wasAccepted  = testr < acceptance; 
  debug_printAccRejc(chain, pf, chain->wasAccepted); 
  chain->prevProposal = pf;   

  if(chain->wasAccepted)
    {
      cntAccept(&(pf->sCtr)); 
      saveAlignAndTreeState(chain); 
    }
  else
    {
      pf->reset_func(chain, pf); 

      cntReject(&(pf->sCtr)); 
      restoreAlignAndTreeState(chain); 

      double diff = fabs(prevLnl - chain->tr->likelihood); 
      if( diff > 1e-6 )
	{
	  printf("problem restoring previous tr/aln state DIFF: %f\n", diff); 
	  assert(0); 
	}
    }

  debug_checkTreeConsistency(chain);
  debug_checkLikelihood(chain); 

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
	      pf->autotune(chain, pf);
	    }
	}
    }
#endif

  chain->currentGeneration++; 
  /* printf("\n");  */
}
