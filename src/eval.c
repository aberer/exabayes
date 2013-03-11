#include "common.h"

#include "axml.h"
/* #include "eval.h"  */

#include "globals.h"

#include "adapterCode.h"

/* call this for verification after the lnl has been evaluated somehow */
static void expensiveVerify(tree *tr)
{  
#ifdef DEBUG_LNL_VERIFY
  double val1 = tr->likelihood; 

  for(int i = 0; i < getNumberOfPartitions(tr) ;++i)
    tr->executeModel[i] = TRUE; 

  evaluateGeneric(tr, tr->start, TRUE); 

  if(processID == 0)
    {
      if(fabs (tr->likelihood - val1 ) > 0.1)
      printf("WARNING: found in expensive evaluation: likelihood difference is %f (with before/after)\t%f\t%f\n", fabs (tr->likelihood - val1 ), val1, tr->likelihood); 
      assert(fabs (tr->likelihood - val1 ) < 0.1);   
    }
  
#endif
}




/** @brief the same as below, but just for one partition 
 */
void evaluateOnePartition(tree *tr, nodeptr start, boolean fullTraversal, int model)
{
  int numPartitions = GET_NUM_PARTITIONS(tr); 

  double perPartitionLH[numPartitions] ; 

  for(int i = 0; i < numPartitions; ++i)
    perPartitionLH[i] = GET_PLH(tr,i); 

  for(int i = 0; i < numPartitions; ++i)
    GET_EXECMODEL(tr, i) = FALSE; 
  GET_EXECMODEL(tr,model) = TRUE; 


#if HAVE_PLL == 1 
  evaluateGeneric(tr, &(gAInfo.partitions), start, fullTraversal); 
#else 
  evaluateGeneric(tr, start, fullTraversal); 
#endif
  
  perPartitionLH[model] = GET_PLH(tr, model);
  for(int i = 0; i < numPartitions; ++i)
    GET_PLH(tr, i) = perPartitionLH[i]; 

  tr->likelihood = 0; 
  for(int i = 0; i < numPartitions; ++i)
    {	
      tr->likelihood +=  GET_PLH(tr,i);
      GET_EXECMODEL(tr,i) = TRUE; 
    }

  expensiveVerify(tr);
}





/** @brief
   only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(tree *tr, nodeptr start, boolean fullTraversal, boolean *models)
{  
  int numPartitions = GET_NUM_PARTITIONS(tr); 
  double perPartitionLH[numPartitions] ; 

  
  for(int i = 0; i < numPartitions; ++i)
    {
      perPartitionLH[i] = GET_PLH(tr,i); 
      GET_EXECMODEL(tr,i) = models[i] ; 
    }

#if HAVE_PLL == 1
  evaluateGeneric(tr, &(gAInfo.partitions), start, fullTraversal); 
#else 
  evaluateGeneric(tr,start, fullTraversal); 
#endif

  /*  correct for hidden examl feature: reduction is applied multiple times */
  for(int i = 0; i < numPartitions; ++i)
    if(models[i] == TRUE)
      perPartitionLH[i] = GET_PLH(tr,i); 

  tr->likelihood = 0; 
  for(int i = 0;i < numPartitions; ++i)
    {
      GET_PLH(tr,i) = perPartitionLH[i]; 
      tr->likelihood += GET_PLH(tr,i); 
      GET_EXECMODEL(tr,i) = TRUE; 
    }

  expensiveVerify(tr);
}


void evaluateGenericWrapper(tree *tr, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL == 1 
  evaluateGeneric(tr, &(gAInfo.partitions), start, fullTraversal);  
#else 
  evaluateGeneric(tr, start, fullTraversal);
#endif

  expensiveVerify(tr);  
}
