
#include "axml.h"
#include "bayes.h"
#include "globals.h"
#include "adapters.h"		

#include "LnlRestorer.hpp"
#include "TreeAln.hpp" 

// #define DEBUG_EVAL


/* call this for verification after the lnl has been evaluated somehow */
static void expensiveVerify(state *chain)
{  
#ifdef DEBUG_LNL_VERIFY
  TreeAln *debugTraln = gAInfo.debugTree;   
  *debugTraln = *(chain->traln); 

  double toVerify = chain->traln->getTr()->likelihood; 
  
  state  *helpChain = (state*) exa_calloc(1,sizeof(chain)) ; 
  helpChain->traln = debugTraln; 
  exa_evaluateGeneric(helpChain, helpChain->traln->getTr()->start, TRUE); 
  double verifiedLnl =  helpChain->traln->getTr()->likelihood; 


  if(chain->currentGeneration != 0 && processID == 0)
    {
      if(fabs (verifiedLnl - toVerify ) > 1e-6)
	printf("WARNING: found in expensive evaluation: likelihood difference is %f (with before/after)\t%f\t%f\n", fabs (verifiedLnl - toVerify ), toVerify, verifiedLnl); 
      assert(fabs (verifiedLnl - toVerify ) < 1e-6);   
    }  

  exa_free(helpChain); 
#endif
}

 
/**
   @brief backs up arrays and executes the newView 
 */ 
void newViewGenericWrapper(state *chain, nodeptr p, boolean masked)
{
  int numberToExecute = 0; 
  int modelToEval = ALL_MODELS; 
  for(int i = 0; i < chain->traln->getNumberOfPartitions(); ++i)
    {
      if(chain->traln->accessExecModel(i))
	{
	  numberToExecute++; 
	  modelToEval = i; 
	}
    }

#ifdef DEBUG_EVAL
  cout << "newViewGenericWrapper on " << p->number  << " and " <<  modelToEval << endl; 
#endif
  
  assert(numberToExecute == 1 || numberToExecute == chain->traln->getNumberOfPartitions()); 
  if(numberToExecute > 1)
    modelToEval = ALL_MODELS; 


  chain->restorer->traverseAndSwitchIfNecessary(p, modelToEval, false); 
  exa_newViewGeneric(chain,p,masked); // NEEDED
}


void evaluateFullNoBackup(state *chain)
{
#ifdef DEBUG_EVAL
  cout << "conducting full evaluation, no backup created" << endl; 
#endif
  
  exa_evaluateGeneric(chain,chain->traln->getTr()->start,TRUE );   
  expensiveVerify(chain);
}




void evaluateGenericWrapper(state *chain, nodeptr start, boolean fullTraversal)
{
#ifdef DEBUG_EVAL
  cout << "evaluateGeneric at " << start->number << "/" << start->back->number << " with " << (fullTraversal ? "TRUE" : "FALSE" ) ; 
#endif

  int model = ALL_MODELS; 
  int numModels = 0; 
  for(int i = 0; i < chain->traln->getNumberOfPartitions() ; ++i)
    {
      if(chain->traln->accessExecModel(i))
	{	  
	  numModels++; 
	  model = i; 
	}
    }
  // cout << "\t" << numModels << "\t"; 
  assert(numModels == 1 || numModels == chain->traln->getNumberOfPartitions()); 
  if(numModels > 1)
    {
      model = ALL_MODELS; 
#ifdef DEBUG_EVAL
      cout << " and all models" << endl; 
#endif
    }
#ifdef DEBUG_EVAL
  else
    cout << " and model "<< model << endl; 
#endif


  chain->restorer->traverseAndSwitchIfNecessary(start, model, fullTraversal);
  chain->restorer->traverseAndSwitchIfNecessary(start->back, model, fullTraversal);

  exa_evaluateGeneric(chain,start,fullTraversal);   
  expensiveVerify(chain);
}


void orientationPointAway(tree *tr, nodeptr p)
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
  
  TreeAln *traln = chain->traln;
  tree *tr = chain->traln->getTr(); 
  int numPartitions = chain->traln->getNumberOfPartitions(); 

  double *perPartitionLH; 
  perPartitionLH = new double[numPartitions]; 

  for(int i = 0; i < numPartitions; ++i)
    perPartitionLH[i] = traln->accessPartitionLH(i); 

  for(int i = 0; i < numPartitions; ++i)
    traln->accessExecModel(i) = FALSE; 
  traln->accessExecModel(model) = TRUE; 

  orientationPointAway(tr, start); 
  orientationPointAway(tr, start->back); 
  
  /* compensating for the fact, that we need to have a tip for full traversal  */
  newViewGenericWrapper(chain, start, TRUE); 

  for(int i = 0; i < numPartitions; ++i)
    traln->accessExecModel(i) = FALSE; 
  traln->accessExecModel(model) = TRUE; 

  newViewGenericWrapper(chain, start->back, TRUE); 

  for(int i = 0; i < numPartitions; ++i)
    traln->accessExecModel(i) = FALSE; 
  traln->accessExecModel(model) = TRUE; 

  evaluateGenericWrapper(chain,start, FALSE);

  perPartitionLH[model] = traln->accessPartitionLH(model) ; 
  for(int i = 0; i < numPartitions; ++i)
    traln->accessPartitionLH(i) = perPartitionLH[i]; 
    

  tr->likelihood = 0; 
  for(int i = 0; i < numPartitions; ++i)
    {	
      tr->likelihood += traln->accessPartitionLH(i);
      traln->accessExecModel(i) = TRUE; 
    }

  delete [] perPartitionLH; 
  expensiveVerify(chain);
}




#if 0 
/**
   @brief only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(state *chain, nodeptr start, boolean fullTraversal, boolean *models)
{  
  tree *tr = chain->traln->getTr(); 
  assert(0);  			/* not in use  */
  
  int numPartitions = getNumberOfPartitions(tr); 
  double *perPartitionLH = new double[numPartitions]; 

  
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


  delete []  perPartitionLH; 
  expensiveVerify(chain); 
}
#endif

