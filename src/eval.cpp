#include <cassert>
#include "axml.h"
#include "GlobalVariables.hpp"
#include "LnlRestorer.hpp"
#include "TreeAln.hpp" 
#include "LikelihoodEvaluator.hpp" 

void evaluateFullNoBackup(TreeAln& traln); 

void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
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

  LikelihoodEvaluator eval; 
  Branch root ; 
  eval.findVirtualRoot(traln, root); 

  nodeptr
    p = root.findNodePtr(debugTraln ); 

  eval.disorientTree(debugTraln, root); 
  exa_evaluateGeneric(debugTraln, p , FALSE); 

  double verifiedLnl =  debugTraln.getTr()->likelihood; 

  if(fabs (verifiedLnl - toVerify ) > ACCEPTED_LIKELIHOOD_EPS)
    {
      tout << "WARNING: found in expensive evaluation: likelihood difference is " 
	   << setprecision(8) <<   fabs (verifiedLnl - toVerify )
	   << " (with toVerify= " << toVerify << ", verified=" << verifiedLnl << ")" << endl; 

      tout << "current tree: " << traln << endl; 
      tout << "help tree: " <<  debugTraln << endl; 	  

      evaluateFullNoBackup(traln);
      tout << "full evaluation on original tree yields: "  << traln.getTr()->likelihood << endl; 
    }  
  assert(fabs (verifiedLnl - toVerify ) < ACCEPTED_LIKELIHOOD_EPS);   
#endif
}


void evaluateFullNoBackup(TreeAln& traln)
{
#ifdef DEBUG_EVAL  
  cout << "conducting full evaluation, no backup created" << endl; 
#endif
  
  exa_evaluateGeneric(traln,traln.getTr()->start,TRUE );   
  expensiveVerify(traln);
}


void exa_newViewParsimony(TreeAln &traln, nodeptr p)
{
#if HAVE_PLL != 0
  newviewParsimony(traln.getTr(), traln.getPartitionsPtr(), p); 
#else 
  newviewParsimony(traln.getTr(),p); 
  // assert(0);
#endif
}


void exa_evaluateParsimony(TreeAln &traln, nodeptr p, boolean fullTraversal, vector<nat> &partitionParsimony)
{
  partitionParsimony.clear();   
  nat* tmp = (nat*)exa_calloc(traln.getNumberOfPartitions(), sizeof(nat)); 
  
#if HAVE_PLL != 0 
  evaluateParsimony(traln.getTr(), traln.getPartitionsPtr(), p, fullTraversal, tmp); 
#else 
  evaluateParsimony(traln.getTr(), p, fullTraversal, tmp); 
  // assert(0);
#endif

  for(int i = 0; i < traln.getNumberOfPartitions(); ++i)
    partitionParsimony.push_back(tmp[i]); 

  exa_free(tmp); 
}


void evaluatePartialNoBackup(TreeAln& traln, nodeptr p)
{  
  exa_evaluateGeneric(traln,p,FALSE );   
  expensiveVerify(traln);
}
