
#include "globals.h"
#include "common.h"

#include "config.h"
#include "axml.h"

#include "proposalStructs.h"
#include "main-common.h"

#include "chain.h"
#include "proposals.h"
#include "output.h"
#include "convergence.h" 



void initDefaultValues(state *theState, tree *tr); 
void addInitParameters(state *curstate, initParamStruct *initParams); 


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



void initializeIndependentChains(tree *tr, state **resultIndiChains, initParamStruct **initParamsPtr)
{
#ifndef _USE_NCL_PARSER
  /* we MUST use the ncl parser currently */
  assert(0); 
#endif

  FILE *treeFH = NULL; 
  if( numberOfStartingTrees > 0 )
    treeFH = myfopen(tree_file, "r"); 

  *initParamsPtr = calloc(1,sizeof(initParamStruct)); 
  /* initParamStruct *initParams = *initParamsPtr;  */
  
  parseConfigWithNcl(configFileName, initParamsPtr);

  printf("num indi chains in init is %d\n", (*initParamsPtr)->numIndiChains); 

  *resultIndiChains = calloc( (*initParamsPtr)->numIndiChains, sizeof(state)); 

  unsigned int bvLength = 0;   
  tr->bitVectors = initBitVector(tr->mxtips, &bvLength);
  hashtable *ht = initHashTable(tr->mxtips * tr->mxtips * 10);
  
  for(int i = 0; i < (*initParamsPtr)->numIndiChains; ++i)
    {
      state *theState = (*resultIndiChains) + i; 
      theState->id = i; 
      initDefaultValues(theState, tr);
      addInitParameters(theState, *initParamsPtr); 
      normalizeProposalWeights(theState); 

      theState->bvHash = ht; 
      theState->tr = tr; 

      /* init the param dump  */
      theState->dump.topo = malloc(sizeof(bestlist)); 
      theState->dump.topo->ninit = 0;       
      initBestTree(theState->dump.topo, 1, tr->mxtips);
      theState->dump.infoPerPart = calloc(tr->NumberOfModels, sizeof(perPartitionInfo)); 
      theState->dump.branchLengths = calloc(2 * tr->mxtips, sizeof(double)); 

      if( i < numberOfStartingTrees )
	{
	  treeReadLen(treeFH, tr, TRUE, TRUE, FALSE); 
	  if(processID == 0)
	    printf("initializing chain %d with provided starting tree\n", i); 
	}
      else 
	{
	  makeRandomTree(tr);
	  if(processID == 0)
	    printf("initializing chain %d with random tree\n", i); 
	  
	  /* TODO maybe apply prior  */

	  int count = 0; 
	  traverseInitFixedBL( tr->start->back, &count, tr, INIT_BRANCH_LENGTHS );
	  assert(count  == 2 * tr->mxtips - 3);
	}

      /* TODO these are dummy values, we can do better */
      tr->start = find_tip(tr->start, tr );

      /* now save the tree to a chain state */
      saveTreeStateToChain(theState, tr); 
      
      if(processID == 0)
	{	  
	  makeChainFileNames(theState, i); 
	  initializeOutputFiles(theState); 
	}
    }  

  if(numberOfStartingTrees > 0)
    fclose(treeFH); 

  if(processID == 0)
    printf("initialized %d independent chains\n", (*initParamsPtr)->numIndiChains ); 
}



void traverseInitFixedBL(nodeptr p, int *count, tree *tr,  double z )
{
  nodeptr q;
  int i;
  
  for( i = 0; i < tr->numBranches; i++)
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


void traverseAndTreatBL(node *p, tree *tr, double *blBuf, int* cnt, boolean restore)
{
  nodeptr q; 
  assert(tr->numBranches == 1); 

  if(restore == TOPO_RESTORE )
    p->z[0] = blBuf[*cnt]; 
  else if(restore == TOPO_SAVE)    
    blBuf[*cnt] =  p->z[0]; 
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
  /* TODO  */
  assert(tr->numBranches == 1); 

  /* restore topology */
  recallBestTree(chain->dump.topo, 1,tr); 

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 
  
  /* restore model parameters */
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i ; 
      tr->partitionData[i].alpha = info->alpha;
      makeGammaCats(tr->partitionData[i].alpha, tr->partitionData[i].gammaRates, 4, tr->useMedian);

      memcpy(tr->partitionData[i].substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(tr->partitionData[i].frequencies, info->frequencies, 4 * sizeof(double)); 

      initReversibleGTR (tr, i); 
    } 

  evaluateGeneric(tr, tr->start, TRUE ); 

  if(processID == 0)
    printf("switching to chain %d, lnl=%f\n", chain->id, tr->likelihood); 
  
  if( chain->currentGeneration != 0 &&  fabs (tr->likelihood - chain->dump.likelihood ) > LIKELIHOOD_EPSILON )
    {
      if(processID == 0)
	printf("WARNING: obtained a different likelihood (eps: %f) after restoring previous chain state (before/after): %f / %f\n", LIKELIHOOD_EPSILON, chain->dump.likelihood, tr->likelihood); 
      assert( fabs(chain->dump.likelihood - tr->likelihood )  > LIKELIHOOD_EPSILON ) ; 
    }
}


void saveTreeStateToChain(state *chain, tree *tr)
{
  chain->dump.likelihood = tr->likelihood;   

#if HAVE_PLL == 1 
  saveBestTree(chain->dump.topo, tr); 
#else 
  saveBestTree(chain->dump.topo, tr, FALSE); 
#endif

  /* save branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_SAVE); 
  assert(cnt == 2 * tr->mxtips - 3 ); 

  /* save model parameters */
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i; 
      info->alpha =  tr->partitionData[i].alpha; 
      
      memcpy(info->substRates, tr->partitionData[i].substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, tr->partitionData[i].frequencies, 4 * sizeof(double)); 
    }  
}
