
#include "common.h"

#include "config.h"
#include "axml.h"

#include "proposalStructs.h"

#include "globals.h"

#include "main-common.h"

#include "chain.h"
#include "proposals.h"
#include "output.h"
#include "convergence.h" 

#include "treeRead.h"

#include "bayes-topo.h"

#include "randomness.h"
#include "eval.h"
#include "adapterCode.h"


/* #define DEBUG_BL */

void initDefaultValues(state *theState, tree *tr); 
void addInitParameters(state *curstate, initParamStruct *initParams); 

void exa_makeRandomTree(tree *tr);  

void printInfo(state *chain, const char *format, ...)
{  
  if(processID == 0)
    {
      printf("[run %d / heat %d / gen %d] ", chain->id / gAInfo.numberCoupledChains, chain->couplingId, chain->currentGeneration); 
      va_list args;
      va_start(args, format);     
      vprintf(format, args );
      va_end(args);
    }
}


void makeChainFileNames(state *theState, int num)
{
  char tName[1024],
    pName[1024] ; 
  
  sprintf(tName, "%s%s_topologies.%s.%d", workdir, PROGRAM_NAME, run_id, num); 
  sprintf(pName, "%s%s_parameters.%s.%d", workdir, PROGRAM_NAME, run_id, num); 
  
  /* todo binary state file?  */
  theState->topologyFile = fopen(tName, "w"); 
  theState->outputParamFile = fopen(pName, "w");   
}


double getChainHeat(state *chain )
{
  const double  deltaT = HEAT_FACTOR; 

  if(chain->couplingId == 0 )
    return 1; 
  
  double tmp  = 1. + deltaT * chain->couplingId; 
  
  assert(tmp > 1); 
  double myHeat = 1. / (double)tmp; 

  assert(myHeat < 1.);
  return myHeat; 
}


/* this function needs to be called when branch lengths are read from
   a file. If this is the case, these BLs are not transformed
   correctly. This function corrects for that.  */
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




static  void copyParamDump(tree *tr,paramDump *dest, const paramDump *src)
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



#if HAVE_PLL == 1 
void initializeTree(tree *tr, partitionList *partitions, analdef *adef); 
#else  
void initializeTree(tree *tr, analdef *adef); 
#endif


void preInitTree(tree *tr)
{
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA; 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->saveMemory = FALSE;
  tr->manyPartitions = FALSE;
  tr->categories             = 25;
  tr->grouped = FALSE;
  tr->constrained = FALSE;
  tr->gapyness               = 0.0; 
  tr->useMedian = FALSE;
}


void initializeIndependentChains(tree *tr, analdef *adef, state **resultIndiChains, initParamStruct **initParamsPtr)
{
  gAInfo.successFullSwitchesBatch = 0; 

  FILE *treeFH = NULL; 
  if( gAInfo.numberOfStartingTrees > 0 )
    treeFH = myfopen(tree_file, "r"); 

  *initParamsPtr = exa_calloc(1,sizeof(initParamStruct)); 
  /* initParamStruct *initParams = *initParamsPtr;  */

  parseConfigWithNcl(configFileName, initParamsPtr);

  gAInfo.numberOfRuns =   (*initParamsPtr)->numIndiChains; 
  gAInfo.numberCoupledChains = (*initParamsPtr)->numCoupledChains; 
  int totalNumChains = gAInfo.numberOfRuns * gAInfo.numberCoupledChains; 

#ifdef MC3_SPACE_FOR_TIME
  int treesNeeded = gAInfo.numberCoupledChains  -1 ; 
#if HAVE_PLL == 1 
  gAInfo.partitions = (partitionList**)exa_realloc(gAInfo.partitions, (treesNeeded + 1)  * sizeof(partitionList*)); 
  for(int i = 1 ; i < treesNeeded+1; ++i)
    gAInfo.partitions[i] = exa_calloc(1,sizeof(partitionList));
#endif
  tree *moreTrees = exa_calloc(treesNeeded, sizeof(tree)); 
  for(int i = 0; i < treesNeeded; ++i)
    {
      preinitTree(moreTrees+i); 
#if HAVE_PLL == 1 
      initializeTree(moreTrees + i, gAInfo.partitions[i+1] ,adef); 
#else 
      initializeTree(moreTrees + i ,adef); 
#endif
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
  
  for(int i = 0; i < totalNumChains; ++i)
    { 
      state *theChain = *resultIndiChains + i;     
      theChain->bvHash = ht; 

      theChain->id = i; 
      theChain->couplingId = i % gAInfo.numberCoupledChains ; 
      
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
      tree *myTree = theChain->tr; 
      myTree->bitVectors = tr->bitVectors; 

      initDefaultValues(theChain, myTree);
      addInitParameters(theChain, *initParamsPtr); 
      normalizeProposalWeights(theChain); 

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

  if( chain->dump.likelihood != 0 && fabs (tr->likelihood - chain->dump.likelihood ) > 0.000001 )
    {
      printInfo(chain, "WARNING: obtained a different likelihood  after restoring previous chain state (before/after): %f / %f\n", 
		chain->id, chain->currentGeneration, chain->dump.likelihood, tr->likelihood); 
      assert( fabs(chain->dump.likelihood - tr->likelihood ) < 0.000001 ) ; 
    }
}



/**
   @brief Save all relevan information from the tree into the chain
   state. 
 */ 
void saveTreeStateToChain(state *chain)
{
  tree *tr  = chain->tr; 
  chain->dump.likelihood = tr->likelihood;   
  
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
