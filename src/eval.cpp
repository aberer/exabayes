#include <cassert>
#include "axml.h"
#include "GlobalVariables.hpp"
#include "output.h"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp" 


static void exa_newViewGeneric(TreeAln& traln, nodeptr p, boolean masked)
{
#if HAVE_PLL != 0
  newviewGeneric(traln.getTr(), traln.getPartitionsPtr(), p, masked); 
#else 
  newviewGeneric(traln.getTr(), p, masked); 
#endif 
}


static void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL != 0
  evaluateGeneric(traln.getTr(), traln.getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(traln.getTr(), start, fullTraversal); 
#endif  
}

// grml 
#if HAVE_PLL == 0
extern "C"
{
  void newviewParsimony(tree *tr, nodeptr  p); 
  nat evaluateParsimony(tree *tr, nodeptr p, boolean full, nat *partitionParsimony); 
}
#else 
extern "C"
{
  void newviewParsimony(tree *tr, partitionList *pr, nodeptr  p); 
  nat evaluateParsimony(tree *tr, partitionList *pr, nodeptr p, boolean full, nat *partitionParsimony); 
}
#endif

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


/* call this for verification after the lnl has been evaluated somehow */
void expensiveVerify(TreeAln& traln)
{  
#ifdef DEBUG_LNL_VERIFY
  TreeAln &debugTraln = *(globals.debugTree);   
  double toVerify = traln.getTr()->likelihood; 

  debugTraln = traln; 
  
#if 1 
  nodeptr
    p = findNodeFromBranch(debugTraln.getTr(), findRoot(traln.getTr())); 
  orientationPointAway(debugTraln.getTr(), p);
  orientationPointAway(debugTraln.getTr(), p->back);

  exa_evaluateGeneric(debugTraln, p , FALSE); 
#else 
  exa_evaluateGeneric(debugTraln, debugTraln.getTr()->start , TRUE); 
#endif
  double verifiedLnl =  debugTraln.getTr()->likelihood; 


  if( isOutputProcess())
    {      
      if(fabs (verifiedLnl - toVerify ) > ACCEPTED_LIKELIHOOD_EPS)
	{
	  tout << "WARNING: found in expensive evaluation: likelihood difference is " 
	       << setprecision(8) <<   fabs (verifiedLnl - toVerify )
	       << "(with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")" << endl; 

	  tout << "current tree: " << traln << endl; 
	  tout << "help tree: " <<  debugTraln << endl; 
	  
	}
      assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS);   
    }  
#endif
}


void exa_newViewParsimony(TreeAln &traln, nodeptr p)
{
#if HAVE_PLL != 0
  newviewParsimony(traln.getTr(), traln.getPartitionsPtr(), p); 
#else 
  // newviewParsimony(traln.getTr(),p); 
  assert(0);
#endif
}



void exa_evaluateParsimony(TreeAln &traln, nodeptr p, boolean fullTraversal, vector<nat> &partitionParsimony)
{
  partitionParsimony.clear();   
  // nat tmp[traln.getNumberOfPartitions()]; 
  nat* tmp = (nat*)exa_calloc(traln.getNumberOfPartitions(), sizeof(nat)); 
  
#if HAVE_PLL != 0 
  evaluateParsimony(traln.getTr(), traln.getPartitionsPtr(), p, fullTraversal, tmp); 
#else 
  // return evaluateParsimony(traln.getTr(),  p, fullTraversal); 
  assert(0);
  // return 0 ;
#endif

  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    partitionParsimony.push_back(tmp[i]); 

  exa_free(tmp); 
}

 
/**
   @brief backs up arrays and executes the newView 
 */ 
void newViewGenericWrapper(TreeAln &traln, nodeptr p, boolean masked)
{
  int numberToExecute = 0; 
  int modelToEval = ALL_MODELS; 
  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    {
      if(traln.accessExecModel(i))
	{
	  numberToExecute++; 
	  modelToEval = i; 
	}
    }

  assert(numberToExecute == 1 || numberToExecute == traln.getNumberOfPartitions()); 
  if(numberToExecute > 1)
    modelToEval = ALL_MODELS; 

#ifdef DEBUG_EVAL
  if(isOutputProcess())
  cout << "newViewGenericWrapper on " << p->number  << " and model " <<  modelToEval << endl; 
#endif

  if(p->x)
    {
      tree *tr = traln.getTr();
      assert(NOT isTip(p->number, tr->mxtips)); 
      p->x = 0; 
      p->next->x = 1; 
    }
  traln.getRestorer()->traverseAndSwitchIfNecessary(traln, p, modelToEval, false); 
  exa_newViewGeneric(traln,p,masked); // NEEDED
}



void evaluatePartialNoBackup(TreeAln& traln, nodeptr p)
{  
  exa_evaluateGeneric(traln,p,FALSE );   
  expensiveVerify(traln);
}

void evaluateFullNoBackup(TreeAln& traln)
{
#ifdef DEBUG_EVAL
  if(isOutputProcess())
  cout << "conducting full evaluation, no backup created" << endl; 
#endif
  
  exa_evaluateGeneric(traln,traln.getTr()->start,TRUE );   
  expensiveVerify(traln);
}




void evaluateGenericWrapper(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
#ifdef DEBUG_EVAL
  if(isOutputProcess())
  cout << "evaluateGeneric at " << start->number << "/" << start->back->number << " with " << (fullTraversal ? "TRUE" : "FALSE" ) ; 
#endif

  int model = ALL_MODELS; 
  int numModels = 0; 
  for(int i = 0; i < traln.getNumberOfPartitions() ; ++i)
    {
      if(traln.accessExecModel(i))
	{	  
	  numModels++; 
	  model = i; 
	}
    }
  assert(numModels == 1 || numModels == traln.getNumberOfPartitions()); 
  if(numModels > 1)
    {
      model = ALL_MODELS; 
#ifdef DEBUG_EVAL
      if(isOutputProcess())
      cout << " and all models" << endl; 
#endif
    }
#ifdef DEBUG_EVAL
  else
    if(isOutputProcess())
    cout << " and model "<< model << endl; 
#endif

  traln.getRestorer()->traverseAndSwitchIfNecessary(traln, start, model, fullTraversal);
  traln.getRestorer()->traverseAndSwitchIfNecessary(traln, start->back, model, fullTraversal);
  
  exa_evaluateGeneric(traln,start,fullTraversal);   
#ifdef DEBUG_LNL_VERIFY
  if(globals.verifyLnl)
    expensiveVerify(traln);
#endif
}




void evaluatePartitions(TreeAln &traln, nodeptr start, vector<nat> models  ) 
{
  // notice assuming full traversal -- need partial for per-gene branch lengths 

  tree *tr = traln.getTr(); 
  nat numPart = traln.getNumberOfPartitions();
  
  vector<double> perPartitionLH; 
  for(auto i = 0; i < numPart; ++i)
    perPartitionLH.push_back(traln.accessPartitionLH(i));
  
  for(auto i = 0; i < numPart; ++i)
    traln.accessExecModel(i) = FALSE; 
  for(auto m : models)
    traln.accessExecModel(m) = TRUE ; 
  

  // SERIOUSLY??? 
#ifdef DEBUG_LNL_VERIFY
  globals.verifyLnl = false; 	// HACK
#endif
  evaluateGenericWrapper(traln, start, TRUE );
#ifdef DEBUG_LNL_VERIFY
  globals.verifyLnl = true; 	// HACK
#endif
  
  for(auto m : models )
    perPartitionLH[m] = traln.accessPartitionLH(m); 

  for(auto i = 0; i < numPart; ++i )
    traln.accessPartitionLH(i) = perPartitionLH[i]; 

  tr->likelihood = 0; 
  for(auto i = 0; i < numPart; ++i)
    {
      tr->likelihood += traln.accessPartitionLH(i); 
      traln.accessExecModel(i) = TRUE; 
    }  
  
  expensiveVerify(traln);   
}




/**
   @brief the same as below, but just for one partition 

   updates chain likelihood and chain-partition-lnl accordingly. 
 */
void evaluateOnePartition(TreeAln& traln, nodeptr start, boolean fullTraversal, int model)
{
  assert(fullTraversal); 	/* partial tarversal does not make sense */

  tree *tr = traln.getTr(); 
  int numPartitions = traln.getNumberOfPartitions(); 

  double *perPartitionLH; 
  perPartitionLH = new double[numPartitions]; 

  for(int i = 0; i < numPartitions; ++i)
    perPartitionLH[i] = traln.accessPartitionLH(i); 

  for(int i = 0; i < numPartitions; ++i)
    traln.accessExecModel(i) = FALSE; 
  traln.accessExecModel(model) = TRUE; 

  orientationPointAway(tr, start); 
  orientationPointAway(tr, start->back); 
  
  /* compensating for the fact, that we need to have a tip for full traversal  */
  newViewGenericWrapper(traln, start, TRUE); 

  for(int i = 0; i < numPartitions; ++i)
    traln.accessExecModel(i) = FALSE; 
  traln.accessExecModel(model) = TRUE; 

  newViewGenericWrapper(traln, start->back, TRUE); 

  for(int i = 0; i < numPartitions; ++i)
    traln.accessExecModel(i) = FALSE; 
  traln.accessExecModel(model) = TRUE; 
#ifdef DEBUG_LNL_VERIFY
  globals.verifyLnl = false; 	// HACK
#endif
  evaluateGenericWrapper(traln,start, FALSE);
#ifdef DEBUG_LNL_VERIFY
  globals.verifyLnl = true  ; 
#endif

  perPartitionLH[model] = traln.accessPartitionLH(model); 
  for(int i = 0; i < numPartitions; ++i)
    traln.accessPartitionLH(i) = perPartitionLH[i]; 

  tr->likelihood = 0; 
  for(int i = 0; i < numPartitions; ++i)
    {	
      tr->likelihood += traln.accessPartitionLH(i);
      traln.accessExecModel(i) = TRUE; 
    }

  delete [] perPartitionLH; 
  expensiveVerify(traln);
}




#if 0 
/**
   @brief only evaluate partition given in the execute mask "models"
*/
void evaluatePartitions(Chain *chain, nodeptr start, boolean fullTraversal, boolean *models)
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
