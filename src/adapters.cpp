
#include "axml.h"

#include "bayes.h"
#include "globals.h"
#include "TreeAln.hpp" 


void exa_newViewGeneric(state *chain, nodeptr p, boolean masked)
{
#if HAVE_PLL == 1 
  newviewGeneric(chain->traln->getTr(), chain->traln->getPartitionsPtr(), p, masked); 
#else 
  newviewGeneric(chain->traln->getTr(), p, masked); 
#endif 
} 


void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q)
{
#if HAVE_PLL == 1 
  hookupDefault(p,q); 
#else
  hookupDefault(p,q,tr->numBranches); 
#endif
}

void exa_evaluateGeneric(state *chain, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL == 1 
  evaluateGeneric(chain->traln->getTr(), chain->traln->getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(chain->traln->getTr(), start, fullTraversal); 
#endif  
}


