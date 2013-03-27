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
  double toVerify = chain->likelihood; 

  for(int i = 0; i < getNumberOfPartitions(tr) ;++i)
    setExecModel(chain,i,TRUE); 
  
  exa_evaluateGeneric(chain,tr->start,TRUE); 

  if(chain->currentGeneration != 0 && processID == 0)
    {
      if(fabs (tr->likelihood - toVerify ) > 0.1)
      printf("WARNING: found in expensive evaluation: likelihood difference is %f (with before/after)\t%f\t%f\n", fabs (tr->likelihood - toVerify ), toVerify, tr->likelihood); 
      assert(fabs (tr->likelihood - toVerify ) < 0.1);   
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
  for(int i = tr->mxtips+1; i < 2 * tr->mxtips; ++i)
    {
      nodeptr p = tr->nodep[i]; 
      if(p->x)
	chain->lnl.orientation[ctr] =  0; 
      else if(p->next->x)
	chain->lnl.orientation[ctr] = 1; 
      else if(p->next->next->x)
	chain->lnl.orientation[ctr] = 2; 
      else 
	assert(0); 
      ctr++;
    }

  chain->lnl.start = chain->startNode; 
}


/**
   @brief loads the saved orientation of the x-vectors
 */
void loadOrientation(state *chain)
{
  tree *tr = chain->tr; 
  int ctr = 0; 
  for(int i = tr->mxtips+1; i < 2 * tr->mxtips; ++i)
    {
      nodeptr p = tr->nodep[i]; 
      if(chain->lnl.orientation[ctr] == 0)
	{
	  p->x = 1; p->next->x = p->next->next->x = 0; 
	}
      else if(chain->lnl.orientation[ctr] == 1)
	{
	  p->x = p->next->next->x = 0; p->next->x = 1; 
	}
      else if (chain->lnl.orientation[ctr] == 2)
	{
	  p->x  = p->next->x = 0; p->next->next->x = 1; 
	}
      ctr++; 
    }

  assert(0); /* TODO implement */
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
    memcpy(chain->lnl.vectorsPerPartition[model][i] , xVector[i], sizeof(double) * partition->width); 
}


/**
   @brief loads the lnl arrays 
 */
void loadArray(state *chain, int model)
{
  pInfo *partition = getPartition(chain, model); 
  double **xVector = getXPtr(chain, model); 
  for(int i= 0; i < chain->tr->mxtips-2; ++i)
    memcpy(xVector[i], chain->lnl.vectorsPerPartition[model][i], sizeof(double) * partition->width); 
}



void evaluateGenericWrapper(state *chain, nodeptr start, boolean fullTraversal)
{
  exa_evaluateGeneric(chain,start,fullTraversal); 

  int numPartitions = getNumberOfPartitions(chain->tr); 
  for(int i = 0; i < numPartitions; ++i)
    chain->lnl.partitionLnl[i] = getPLH( chain, i );

  chain->lnl.likelihood = chain->tr->likelihood; 

  expensiveVerify(chain);
}




static void updatePartitionLnl(state *chain, double lnl, int model)
{
  chain->lnl.likelihood -= chain->lnl.partitionLnl[model]; 
  chain->lnl.partitionLnl[model] = lnl; 
  chain->lnl.likelihood += chain->lnl.partitionLnl[model];

  double tmp = 0;
  for(int i = 0; i < getNumberOfPartitions(chain->tr); ++i)
    tmp += chain->lnl.partitionLnl[i]; 
  assert(fabs(tmp - chain->lnl.likelihood) < 1e-6); 
}


/**
   @brief the same as below, but just for one partition 

   updates chain likelihood and chain-partition-lnl accordingly. 
 */
void evaluateOnePartition(state *chain, nodeptr start, boolean fullTraversal, int model)
{
  tree *tr = chain->tr; 
  int numPartitions = getNumberOfPartitions(chain->tr); 

  double perPartitionLH[numPartitions] ; 

  for(int i = 0; i < numPartitions; ++i)
    perPartitionLH[i] = getPLH(chain,i); 

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

  updatePartitionLnl(chain, getPLH(chain, model ), model);  

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

