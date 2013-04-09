
#include "axml.h"
#include "bayes.h"
#include "randomness.h"

/* TODO ask alexis, if we can do this w/o evaluating the entire
   tree */




void multiplyBranchLengthsRecursively(tree *tr, nodeptr p, double multiplier)
{  
  double newZ = pow(p->z[0],multiplier); 
  p->z[0] = p->back->z[0] = newZ; 

  if(isTip(p->number, tr->mxtips))
    return; 

  nodeptr q = p->next; 
  while(p != q)
    {
      multiplyBranchLengthsRecursively(tr, q->back,multiplier); 
      q = q->next; 
    }
}


void applyTLMult(state *chain,  proposalFunction *pf )
{
  double multiplier = drawMultiplier(chain, pf->parameters.multiplier);
  pf->remembrance.multiplier = multiplier; 
  multiplyBranchLengthsRecursively(chain->tr, chain->tr->start->back, multiplier); 
}


void resetTLMult(state *chain, proposalFunction *pf)
{
  multiplyBranchLengthsRecursively(chain->tr, chain->tr->start->back, 1/pf->remembrance.multiplier); 
}



