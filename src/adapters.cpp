
#include "axml.h"

#include "bayes.h"
#include "globals.h"
#include "TreeAln.hpp" 
#include "Chain.hpp" 



void exa_newViewGeneric(Chain *chain, nodeptr p, boolean masked)
{
#if HAVE_PLL != 0
  newviewGeneric(chain->traln->getTr(), chain->traln->getPartitionsPtr(), p, masked); 
#else 
  newviewGeneric(chain->traln->getTr(), p, masked); 
#endif 
} 


void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL != 0
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}

void exa_evaluateGeneric(Chain *chain, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL != 0
  evaluateGeneric(chain->traln->getTr(), chain->traln->getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(chain->traln->getTr(), start, fullTraversal); 
#endif  
}




bool isOutputProcess()
{
#if HAVE_PLL != 0
  return true; 
#else 
  return processID == 0; 
#endif 
  
}
