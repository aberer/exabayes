
#include "common.h"

#include "config.h"
#include "axml.h"

#include "globals.h"

#include "proposalStructs.h"
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



void printInfo(state *chain, const char *format, ...)
{
  if(processID == 0)
    {
      printf("[run %d / heat %d / gen %d] ", chain->id / gAInfo.numberOfRuns, chain->couplingId, chain->currentGeneration); 
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


static node *find_tip( node *n, tree *tr ) 
{
  if( isTip(n->number, tr->mxtips) ) {
    return n;
  } else {
    return find_tip( n->back, tr );
  }  
}



double getChainHeat(state *chain )
{
  return (double)(1. / (1. + HEAT_FACTOR * (double)chain->couplingId)); 
}


/* this function needs to be called when branch lengths are read from
   a file. If this is the case, these BLs are not transformed
   correctly. This function corrects for that.  */
void traverseInitCorrect(nodeptr p, int *count, tree *tr )
{
  nodeptr q;
  int i;
  
  assert(NOT HAS_PERGENE_BL(tr)); 
  
  for( i = 0; i < GET_NUM_BRANCHES(tr); i++)
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


void initParamDump(tree *tr, paramDump *dmp)
{  
  dmp->topo = setupTopol(tr->mxtips); 
  dmp->infoPerPart = calloc(GET_NUM_PARTITIONS(tr), sizeof(perPartitionInfo));
  dmp->branchLengths = calloc(2 * tr->mxtips, sizeof(double));
}


void initializeIndependentChains(tree *tr, state **resultIndiChains, initParamStruct **initParamsPtr)
{
#ifndef _USE_NCL_PARSER
  /* we MUST use the ncl parser currently */
  assert(0); 
#endif


  FILE *treeFH = NULL; 
  if( gAInfo.numberOfStartingTrees > 0 )
    treeFH = myfopen(tree_file, "r"); 

  *initParamsPtr = calloc(1,sizeof(initParamStruct)); 
  /* initParamStruct *initParams = *initParamsPtr;  */
  
  parseConfigWithNcl(configFileName, initParamsPtr);
  
  gAInfo.numberOfRuns =   (*initParamsPtr)->numIndiChains; 
  gAInfo.numberCoupledChains = (*initParamsPtr)->numCoupledChains; 
  int totalNumChains = gAInfo.numberOfRuns * gAInfo.numberCoupledChains; 

  printf("number of independent runs=%d, number of coupled chains per run=%d => total of %d chains \n", gAInfo.numberOfRuns, gAInfo.numberCoupledChains, totalNumChains ); 
  *resultIndiChains = calloc( totalNumChains , sizeof(state));   

  unsigned int bvLength = 0;   
  tr->bitVectors = initBitVector(tr->mxtips, &bvLength);
  hashtable *ht = initHashTable(tr->mxtips * tr->mxtips * 10);

  
  for(int i = 0; i < totalNumChains; ++i)
    { 
      printf("setting up chain %d\n" ,i); 

      state *theChain = *resultIndiChains +   i; 
      theChain->bvHash = ht; 
      theChain->tr = tr; 
      theChain->couplingId = i % (*initParamsPtr)->numCoupledChains ; 
      
      theChain->id = i; 

      initDefaultValues(theChain, tr);
      addInitParameters(theChain, *initParamsPtr); 
      normalizeProposalWeights(theChain); 

      /* init the param dump  */
      initParamDump(tr, &(theChain->dump)); 
      
      for(int i = 0; i < GET_NUM_PARTITIONS(tr); ++i ) 
	exa_initReversibleGTR(tr,i);

      /* init rng */
      theChain->rCtr.v[0] = 0; 
      theChain->rCtr.v[1] = 0; 
      randCtr_t r = drawGlobalRandInt();
      theChain->rKey.v[0] = r.v[0]; 
      theChain->rKey.v[1] = r.v[1]; 
      if(processID == 0)
	printf("initialized chain %d with seed %d,%d\n", theChain->id, theChain->rKey.v[0], theChain->rKey.v[1]); 
      
      
      if( i % gAInfo.numberCoupledChains == 0)
	{
	  /* initialize with new tree  */
	  if( i  / gAInfo.numberCoupledChains < gAInfo.numberOfStartingTrees )
	    {
	      boolean hasBranchLength = readTreeWithOrWithoutBL(tr, treeFH); 

	      if(hasBranchLength)
		{
		  int count = 0; 
		  traverseInitCorrect(tr->start->back, &count, tr ) ;
		  assert(count == 2 * tr->mxtips -3); 
		}
	      else 
		{
		  /* set some standard branch lengths */
		  /* TODO based on prior?   */
		  int count = 0; 
		  traverseInitFixedBL( tr->start->back, &count, tr, INIT_BRANCH_LENGTHS );
		  assert(count  == 2 * tr->mxtips - 3);	  	      
		}

	      if(processID == 0)
		printf("initializing chain %d with provided starting tree\n", i); 
	    }
	  else 
	    {
	      makeRandomTree(tr);
	      if(processID == 0)
		printf("initializing chain %d with random tree\n", i); 

	      /* TODO maybe prior for initial branch lengths   */
	      int count = 0; 
	      traverseInitFixedBL( tr->start->back, &count, tr, INIT_BRANCH_LENGTHS );
	      assert(count  == 2 * tr->mxtips - 3);	  
	    }
	}
      else 
	{
	  /* initialize with previous state */
	  applyChainStateToTree( (*resultIndiChains) + i-1, tr); 

	  if(processID == 0)
	    printf("initialize heated chain %d from previous state\n", i); 
	}

      tr->start = find_tip(tr->start, tr );

      evaluateGenericWrapper(tr, tr->start, TRUE); 

      /* now save the tree to a chain chain */
      saveTreeStateToChain(theChain, tr); 

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
  
  for( i = 0; i < GET_NUM_BRANCHES(tr); i++)
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
  assert(GET_NUM_BRANCHES(tr) == 1); 

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





void applyChainStateToTree(state *chain, tree *tr)
{
  /* TODO enable multi-branch    */
  assert(GET_NUM_BRANCHES(tr) == 1); 

  boolean treeWasRestored = restoreTree(chain->dump.topo, tr); 
  assert(treeWasRestored);   

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 

  /* restore model parameters */
  for(int i = 0; i < GET_NUM_PARTITIONS(tr); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i ; 

      pInfo *partition = & GET_PARTITION(tr, i); 
      partition->alpha = info->alpha;

      memcpy(partition->substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(partition->frequencies, info->frequencies, 4 * sizeof(double));

      exa_initReversibleGTR (tr, i); 

      makeGammaCats(partition->alpha, partition->gammaRates, 4, tr->useMedian);
    } 

  evaluateGenericWrapper(tr, tr->start, TRUE ); 

  /* printInfo(chain, "switching to run %d / heat %d, lnl=%f\n", chain->id / numberOfRuns, chain->couplingId, tr->likelihood);  */
  if(processID == 0)
    {
#ifdef DEBUG_BL
      printf("branch lengths: ");
      int count = 0; 
      traverseAndPrint(tr->start->back, &count, tr); 
      printf("\n"); 
#endif
    }

  if( fabs (tr->likelihood - chain->dump.likelihood ) > 0.000001 )
    {
      printInfo(chain, "WARNING: obtained a different likelihood  after restoring previous chain state (before/after): %f / %f\n", 
		chain->id, chain->currentGeneration, chain->dump.likelihood, tr->likelihood); 
      assert( fabs(chain->dump.likelihood - tr->likelihood ) < 0.000001 ) ; 
    }
}


void saveTreeStateToChain(state *chain, tree *tr)
{
  chain->dump.likelihood = tr->likelihood;   
  
  saveTree(tr, chain->dump.topo);
  /* assert(treeWasSaved);  */

  /* save branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_SAVE); 
  assert(cnt == 2 * tr->mxtips - 3 ); 

  /* save model parameters */
  for(int i = 0; i < GET_NUM_PARTITIONS(tr); ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i; 
      info->alpha =  GET_PARTITION(tr,i).alpha; 
      
      memcpy(info->substRates, GET_PARTITION(tr,i).substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, GET_PARTITION(tr,i).frequencies, 4 * sizeof(double)); 
    }  


  if(processID == 0)
    {

#ifdef DEBUG_BL
      int count = 0;       
      printf("saving state: ");
      traverseAndPrint(tr->start->back,  &count,tr);
      printf("\n");
#endif



      /*     printf("\n\nsaving state: "); */

      /*     for(int i = 0; i < getNumberOfPartitions(tr); ++i) */
      /*     printf("\nalpha: %f\n", info->alpha) */

      /*     printf("\nRATES:") */
      
      /*     for(int i = 0; i < 6; ++i) */
      /* 	printf("%f,", info->substRates[i]); */
      /*     printf("\n"); */

      /*     printf("BL: "); */
      /*     for(int i = 0; i < 2 * tr->mxtips -3 ;++i) */
      /* 	printf("%f,", info->) */

    }
}
