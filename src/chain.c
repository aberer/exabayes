
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

#include "treeRead.h"

#include "bayes-topo.h"



/* #define DEBUG_BL */

void initDefaultValues(state *theState, tree *tr); 
void addInitParameters(state *curstate, initParamStruct *initParams); 



void printInfo(state *chain, const char *format, ...)
{
  if(processID == 0)
    {
      printf("[chain %d / gen %d] ", chain->id, chain->currentGeneration); 
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



/* this function needs to be called when branch lengths are read from
   a file. If this is the case, these BLs are not transformed
   correctly. This function corrects for that.  */
void traverseInitCorrect(nodeptr p, int *count, tree *tr )
{
  nodeptr q;
  int i;
  
  assert(tr->numBranches > 0); 

  for( i = 0; i < tr->numBranches; i++)
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
      theState->dump.topo = setupTopol(tr->mxtips); 

      theState->dump.infoPerPart = calloc(tr->NumberOfModels, sizeof(perPartitionInfo));
      theState->dump.branchLengths = calloc(2 * tr->mxtips, sizeof(double));
      theState->dump.fracchanges = calloc(tr->NumberOfModels, sizeof(double)); 
      
      for(int i = 0; i < tr->NumberOfModels; ++i )	
	initReversibleGTR(tr,i);

      memcpy(theState->dump.fracchanges, tr->fracchanges, sizeof(double) * tr->NumberOfModels); 
      theState->dump.fracchange =  tr->fracchange; 

      if( i < numberOfStartingTrees )
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

      /* TODO these are dummy values, we can do better */
      tr->start = find_tip(tr->start, tr );

      evaluateGeneric(tr, tr->start, TRUE); 

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
  assert(tr->numBranches == 1); 

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
  assert(tr->numBranches == 1); 

  boolean treeWasRestored = restoreTree(chain->dump.topo, tr); 
  assert(treeWasRestored);   

  /* restore branch lengths */
  int cnt = 0; 
  traverseAndTreatBL(tr->start->back, tr, chain->dump.branchLengths, &cnt, TOPO_RESTORE); 
  assert(cnt == 2 * tr->mxtips -3); 

  /* TODO we never verified, if it is necessary to copy the fracchange (just an assumption)  */
  memcpy(tr->fracchanges, chain->dump.fracchanges, sizeof(double) * tr->NumberOfModels);
  tr->fracchange = chain->dump.fracchange;
  
  /* restore model parameters */
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i ; 
      tr->partitionData[i].alpha = info->alpha;

      memcpy(tr->partitionData[i].substRates, info->substRates , 6 * sizeof(double)); 
      memcpy(tr->partitionData[i].frequencies, info->frequencies, 4 * sizeof(double)); 

      initReversibleGTR (tr, i); 

      makeGammaCats(tr->partitionData[i].alpha, tr->partitionData[i].gammaRates, 4, tr->useMedian);
    } 

  evaluateGeneric(tr, tr->start, TRUE ); 

  printInfo(chain, "switching to chain %d, lnl=%f\n", chain->id, tr->likelihood); 
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
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      perPartitionInfo *info = chain->dump.infoPerPart + i; 
      info->alpha =  tr->partitionData[i].alpha; 
      
      memcpy(info->substRates, tr->partitionData[i].substRates, 6 * sizeof(double)); 
      memcpy(info->frequencies, tr->partitionData[i].frequencies, 4 * sizeof(double)); 
    }  

  memcpy(chain->dump.fracchanges, tr->fracchanges, tr->NumberOfModels * sizeof(double));
  chain->dump.fracchange = tr->fracchange;

  if(processID == 0)
    {

#ifdef DEBUG_BL
      int count = 0;       
      printf("saving state: ");
      traverseAndPrint(tr->start->back,  &count,tr);
      printf("\n");
#endif



      /*     printf("\n\nsaving state: "); */

      /*     for(int i = 0; i < tr->NumberOfModels; ++i) */
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
