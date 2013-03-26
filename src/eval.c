/**
   @file eval.h
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








/** @brief
   only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(state *chain, nodeptr start, boolean fullTraversal, boolean *models)
{  
  tree *tr = chain->tr; 
  assert(0); 
  
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

