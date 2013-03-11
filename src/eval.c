#include "common.h"

#include "axml.h"
/* #include "eval.h"  */

#include "globals.h"

/* call this for verification after the lnl has been evaluated somehow */
static void expensiveVerify(tree *tr)
{  
#ifdef DEBUG_LNL_VERIFY
  double val1 = tr->likelihood; 

  for(int i = 0; i < tr->NumberOfModels ;++i)
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
  /* this most likely works differently in the PLL */
  assert(HAVE_PLL == 0); 

  double perPartitionLH[tr->NumberOfModels] ; 
  memcpy(perPartitionLH, tr->perPartitionLH, sizeof(double) * tr->NumberOfModels); 

  memset(tr->executeModel, 0, sizeof(boolean) * tr->NumberOfModels); 
  tr->executeModel[model] = TRUE ; 

  evaluateGeneric(tr, start, fullTraversal); 
  
  perPartitionLH[model] = tr->perPartitionLH[model];
  memcpy(tr->perPartitionLH, perPartitionLH, sizeof(double) * tr->NumberOfModels); 

  tr->likelihood = 0; 
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {	
      tr->likelihood +=  tr->perPartitionLH[i];     
      tr->executeModel[i] = TRUE; 
    }

  expensiveVerify(tr);
}





/** @brief
   only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(tree *tr, nodeptr start, boolean fullTraversal, boolean *models)
{  
  /* this most likely works differently in the PLL */
  assert(HAVE_PLL == 0); 

  double perPartitionLH[tr->NumberOfModels] ; 

  memcpy(perPartitionLH, tr->perPartitionLH, sizeof(double) * tr->NumberOfModels); 
  memcpy(tr->executeModel, models, sizeof(boolean)* tr->NumberOfModels);   

  evaluateGeneric(tr,start, fullTraversal); 

  /*  correct for hidden examl feature: reduction is applied multiple times */
  for(int i = 0; i < tr->NumberOfModels; ++i)
    if(models[i] == TRUE)
      perPartitionLH[i] = tr->perPartitionLH[i]; 
  memcpy(tr->perPartitionLH, perPartitionLH, sizeof(double) * tr->NumberOfModels); 
  
  tr->likelihood = 0; 
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      if(models[i] == FALSE) 
	tr->perPartitionLH[i] = perPartitionLH[i]; 
      tr->likelihood +=  tr->perPartitionLH[i]; 
      tr->executeModel[i] = TRUE; 
    }

  expensiveVerify(tr);
}





void evaluateGenericWrapper(tree *tr, nodeptr start, boolean fullTraversal)
{
  evaluateGeneric(tr, start, fullTraversal);
  expensiveVerify(tr);  
}
