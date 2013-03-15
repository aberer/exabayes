#include "common.h"
#include "axml.h"
#include "proposalStructs.h"

#include "globals.h"

#include "adapterCode.h"

/* call this for verification after the lnl has been evaluated somehow */
static void expensiveVerify(tree *tr)
{  
#ifdef DEBUG_LNL_VERIFY
  double val1 = tr->likelihood; 

  for(int i = 0; i < getNumberOfPartitions(tr) ;++i)
    setExecModel(tr,i,TRUE); 
  
  exa_evaluateGeneric(tr,tr->start,TRUE); 

  if(processID == 0)
    {
      if(fabs (tr->likelihood - val1 ) > 0.1)
      printf("WARNING: found in expensive evaluation: likelihood difference is %f (with before/after)\t%f\t%f\n", fabs (tr->likelihood - val1 ), val1, tr->likelihood); 
      assert(fabs (tr->likelihood - val1 ) < 0.1);   
    }
  
#endif
}



void evaluateGenericWrapper(tree *tr, nodeptr start, boolean fullTraversal)
{
  exa_evaluateGeneric(tr,start,fullTraversal); 

  expensiveVerify(tr);  
}


/** @brief the same as below, but just for one partition 
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

  exa_evaluateGeneric(tr, start, fullTraversal); 
  
  perPartitionLH[model] = getPLH(chain, model);
  for(int i = 0; i < numPartitions; ++i)
    setPLH(chain,i,perPartitionLH[i]); 

  tr->likelihood = 0; 
  for(int i = 0; i < numPartitions; ++i)
    {	
      tr->likelihood +=  getPLH(chain,i);
      setExecModel(chain,i,TRUE); 
    }

  expensiveVerify(tr);
}





/** @brief
   only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(state *chain, nodeptr start, boolean fullTraversal, boolean *models)
{  
  tree *tr = chain->tr; 
  
  int numPartitions = getNumberOfPartitions(tr); 
  double perPartitionLH[numPartitions] ; 

  
  for(int i = 0; i < numPartitions; ++i)
    {
      perPartitionLH[i] = getPLH(chain,i); 
      setExecModel(chain,i,models[i]); 
    }

  exa_evaluateGeneric(tr, start, fullTraversal); 

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

  expensiveVerify(tr);
}

