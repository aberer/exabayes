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
  
  /* assert(0); /\* TODO implement *\/ */
  /* TODO */
  /* tr->start = chain->lnl.start;  */
}



/**
   @brief saves the lnl arrays 
 */ 
void saveArray(state *chain, int model)
{
  pInfo *partition = getPartition(chain, model );
  double **xVector = getXPtr(chain,model); 
  for(int i = 0; i < chain->tr->mxtips-2; ++i)
    memcpy(chain->lnl.vectorsPerPartition[model][i] , xVector[i], sizeof(double) * (partition->upper - partition->lower) * LENGTH_LNL_ARRAY);   
}


/**
   @brief loads the lnl arrays 
 */
void loadArray(state *chain, int model)
{
  pInfo *partition = getPartition(chain, model); 
  double **xVector = getXPtr(chain, model); 
  for(int i= 0; i < chain->tr->mxtips-2; ++i)
    memcpy(xVector[i], chain->lnl.vectorsPerPartition[model][i], sizeof(double) * (partition->upper - partition->lower) * LENGTH_LNL_ARRAY); 
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

  /* printf("state after restore: \n");  */
  /* printAlnTrState(chain);  */

  branch root = findRoot(chain); 

  /* check, if everything is okay  */
  /* printf("RESTORE: root is {%d,%d}\n", root.thisNode, root.thatNode);  */

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


/* static void updatePartitionLnl(state *chain, double lnl, int model) */
/* { */
/*   chain->lnl.likelihood -= chain->lnl.partitionLnl[model];  */
/*   chain->lnl.partitionLnl[model] = lnl;  */
/*   chain->lnl.likelihood += chain->lnl.partitionLnl[model]; */

/*   double tmp = 0; */
/*   for(int i = 0; i < getNumberOfPartitions(chain->tr); ++i) */
/*     tmp += chain->lnl.partitionLnl[i];  */
/*   assert(fabs(tmp - chain->lnl.likelihood) < 1e-6);  */
/* } */



/**
   @brief the same as below, but just for one partition 

   updates chain likelihood and chain-partition-lnl accordingly. 
 */
void evaluateOnePartition(state *chain, nodeptr start, boolean fullTraversal, int model)
{
  tree *tr = chain->tr; 
  int numPartitions = getNumberOfPartitions(chain->tr); 

  double perPartitionLH[numPartitions] ; 
  

  printf("plh before eval %d: ", model); 
  for(int i = 0; i < numPartitions; ++i)
    {
      perPartitionLH[i] = getPLH(chain,i); 
      printf("%f," , perPartitionLH[i]); 
    }
  printf("\n"); 
  
  

  for(int i = 0; i < numPartitions; ++i)
    setExecModel(chain,i,FALSE); 
  setExecModel(chain,model,TRUE); 

  exa_evaluateGeneric(chain, start, fullTraversal); 
  
  perPartitionLH[model] = getPLH(chain, model);
  for(int i = 0; i < numPartitions; ++i)
    setPLH(chain,i,perPartitionLH[i]); 

  tr->likelihood = 0; 
  for(int i = 0; i < numPartitions; ++i)
    {	
      tr->likelihood += getPLH(chain,i);
      setExecModel(chain,i,TRUE); 
    }
  
  printf("plh after eval :"); 
  for(int i = 0; i < numPartitions; ++i)
    printf("%f,", perPartitionLH[i]); 
  printf("\n"); 

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
  branch root = findRoot(chain); 
  printf("root: {%d,%d}\n", root.thisNode,root.thatNode); 
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"); 
}


      
