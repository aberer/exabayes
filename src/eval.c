/**
   @file eval.c
   @brief Functions for likelihood evaluations.  

   New convention: all these functions are responsible for updating
   the chain->likelihood + chain->partitionLnl accordingly. 
*/ 


#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "adapters.h"		

/* call this for verification after the lnl has been evaluated somehow */
static void expensiveVerify(state *chain)
{  
#ifdef DEBUG_LNL_VERIFY
  tree *tr = chain->tr; 
  double toVerify = chain->tr->likelihood; 
  
  for(int i = 0; i < getNumberOfPartitions(tr) ;++i)
    setExecModel(chain,i,TRUE); 
  
  exa_evaluateGeneric(chain,tr->start,TRUE); 
  
  if(chain->currentGeneration != 0 && processID == 0)
    {
      if(fabs (tr->likelihood - toVerify ) > 1e-6)
      printf("WARNING: found in expensive evaluation: likelihood difference is %f (with before/after)\t%f\t%f\n", fabs (tr->likelihood - toVerify ), toVerify, tr->likelihood); 
      assert(fabs (tr->likelihood - toVerify ) < 1e-6);   
    }  
#endif
}



/**
   @brief saves the current orientation of the x-vectors
 */
void saveOrientation(state *chain)
{
  tree *tr = chain->tr; 
  int ctr = 0; 

  for(int i = tr->mxtips+1; i < 2 * tr->mxtips-1; ++i)
    {
      nodeptr p = tr->nodep[i]; 

      int *val = chain->lnl.orientation + ctr; 
      if(p->x)
	*val = p->back->number; 
      else if(p->next->x)
	*val = p->next->back->number; 
      else if(p->next->next->x)
	*val = p->next->next->back->number;       
      ctr ++; 
    }
}


/**
   @brief loads the saved orientation of the x-vectors
 */
void loadOrientation(state *chain)
{
  tree *tr = chain->tr; 
  int ctr = 0; 
  for(int i = tr->mxtips+1; i < 2 * tr->mxtips-1; ++i)
    {
      int val = chain->lnl.orientation[ctr]; 
      branch b = constructBranch(tr->nodep[i]->number, val ); 
      branchExists(chain->tr, b); 
      
      nodeptr q = findNodeFromBranch(chain->tr, b);
      assert(q->back->number == val); 
      q->x = 1; q->next->x = 0; q->next->next->x = 0;
      
      ctr++; 
    }  
}

/**
   @brief saves the lnl arrays 
 */ 
void saveArray(state *chain, int model)
{
  pInfo *partition = getPartition(chain, model );
  double **xVector = getXPtr(chain,model); 
  for(int i = 0; i < chain->tr->mxtips-2; ++i)
    {
#if HAVE_PLL == 1
      int length = (partition->upper - partition->lower); 
#else 
      int length =   partition->width; 
#endif
      memcpy(chain->lnl.vectorsPerPartition[model][i] , xVector[i], sizeof(double) * length * LENGTH_LNL_ARRAY);   
    }
  
  memcpy(chain->lnl.partitionScaler[model], partition->globalScaler, sizeof(nat) * 2 * chain->tr->mxtips); 
}


/**
   @brief loads the lnl arrays 
 */
void loadArray(state *chain, int model)
{
  pInfo *partition = getPartition(chain, model); 
  double **xVector = getXPtr(chain, model); 
  for(int i= 0; i < chain->tr->mxtips-2; ++i)
    {      
#if HAVE_PLL == 1
      int length = (partition->upper - partition->lower); 
#else 
      int length = partition->width; 
#endif
      memcpy(xVector[i], chain->lnl.vectorsPerPartition[model][i], sizeof(double) *  length  * LENGTH_LNL_ARRAY); 
    }
  memcpy(partition->globalScaler, chain->lnl.partitionScaler[model], sizeof(nat) * 2 * chain->tr->mxtips); 
}


/**
   @brief resets the tree orientation and likelihood arrays
 */ 
void restoreAlignAndTreeState(state *chain)
{  
  int numPart = getNumberOfPartitions( chain->tr); 
  loadOrientation(chain);
  for(int i = 0; i < numPart; ++i)
    loadArray(chain,i); 

  branch root = findRoot(chain->tr); 

  nodeptr p = findNodeFromBranch(chain->tr,root);
  evaluateGenericWrapper(chain, p,FALSE);
}


/**
   @brief saves tree orientation and likelihood arrays 
 */ 
void saveAlignAndTreeState(state *chain)
{
  int numPart = getNumberOfPartitions(chain->tr);
  saveOrientation(chain);
  for(int i = 0; i < numPart; ++i)
    saveArray(chain,i); 
}


void evaluateGenericWrapper(state *chain, nodeptr start, boolean fullTraversal)
{
  exa_evaluateGeneric(chain,start,fullTraversal); 
  expensiveVerify(chain);
}




static void orientationPointAway(tree *tr, nodeptr p)
{
  if(NOT isTip(p->number, tr->mxtips))
    {
      if(p->x)
	{
	  p->next->x = 1; 
	  p->x = 0; 
	}
      
      orientationPointAway(tr, p->next->back); 
      orientationPointAway(tr, p->next->next->back); 
    }
}




/**
   @brief the same as below, but just for one partition 

   updates chain likelihood and chain-partition-lnl accordingly. 
 */
void evaluateOnePartition(state *chain, nodeptr start, boolean fullTraversal, int model)
{
  assert(fullTraversal); 	/* partial tarversal does not make sense */

  tree *tr = chain->tr; 
  int numPartitions = getNumberOfPartitions(chain->tr); 

  double perPartitionLH[numPartitions] ; 

  for(int i = 0; i < numPartitions; ++i)
      perPartitionLH[i] = getPLH(chain,i); 

  for(int i = 0; i < numPartitions; ++i)
    setExecModel(chain,i,FALSE); 
  setExecModel(chain,model,TRUE); 

  orientationPointAway(tr, start); 
  orientationPointAway(tr, start->back); 
  
  /* compensating for the fact, that we need to have a tip for full traversal  */
  exa_newViewGeneric(chain, start, TRUE); 
  exa_newViewGeneric(chain, start->back, TRUE); 
  evaluateGenericWrapper(chain,start, FALSE);

  perPartitionLH[model] = getPLH(chain, model);
  for(int i = 0; i < numPartitions; ++i)
    setPLH(chain,i,perPartitionLH[i]); 

  tr->likelihood = 0; 
  for(int i = 0; i < numPartitions; ++i)
    {	
      tr->likelihood += getPLH(chain,i);
      setExecModel(chain,i,TRUE); 
    }

  expensiveVerify(chain);
}



/**
   @brief only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(state *chain, nodeptr start, boolean fullTraversal, boolean *models)
{  
  tree *tr = chain->tr; 
  assert(0);  			/* not in use  */
  
  int numPartitions = getNumberOfPartitions(tr); 
  double perPartitionLH[numPartitions] ; 

  
  for(int i = 0; i < numPartitions; ++i)
    {
      perPartitionLH[i] = getPLH(chain,i); 
      setExecModel(chain,i,models[i]); 
    }

  exa_evaluateGeneric(chain, start, fullTraversal); 

  /*  correct for hidden examl feature: reduction is applied multiple times */
  for(int i = 0; i < numPartitions; ++i)
    if(models[i] == TRUE)
      perPartitionLH[i] = getPLH(chain,i); 

  tr->likelihood = 0; 
  for(int i = 0;i < numPartitions; ++i)
    {
      setPLH(chain,i,perPartitionLH[i]); 
      tr->likelihood += getPLH(chain,i); 
      setExecModel(chain,i,TRUE); 
    }

  expensiveVerify(chain);      

}



void printAlnTrState(state *chain)
{
  tree *tr = chain->tr; 
  int numPart = getNumberOfPartitions(tr); 
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\norient: "); 
  for(int i = 0; i < tr->mxtips-2; ++i)
    printf("%d,", chain->lnl.orientation[i]); 
  printf("\nlnlvect: \n");
  for(int i = 0; i < numPart; ++i)
    {
      pInfo *partition = getPartition(chain,i);
      int length = (partition->upper-partition->lower) * LENGTH_LNL_ARRAY;

      for(int j = 0; j < tr->mxtips-2; ++j)
  	{
  	  printf("%d: ", j);
  	  for(int k = 0; k < length; ++k)
  	    printf("%f,", chain->lnl.vectorsPerPartition[i][j][k]);
  	  printf("\n");
  	}
    }
  branch root = findRoot(chain->tr); 
  printf("root: {%d,%d}\n", root.thisNode,root.thatNode); 
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"); 
}

