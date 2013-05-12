#include "axml.h"

// #include "proposalFunction.h"
#include "GlobalVariables.hpp"
#include "TreeAln.hpp" 

void exa_newViewGeneric(TreeAln& traln, nodeptr p, boolean masked)
{
#if HAVE_PLL != 0
  newviewGeneric(traln.getTr(), traln.getPartitionsPtr(), p, masked); 
#else 
  newviewGeneric(traln.getTr(), p, masked); 
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



void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal)
{
#if HAVE_PLL != 0
  evaluateGeneric(traln.getTr(), traln.getPartitionsPtr(), start, fullTraversal); 
#else 
  evaluateGeneric(traln.getTr(), start, fullTraversal); 
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
