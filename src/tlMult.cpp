
#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "TreeAln.hpp"

/* TODO ask alexis, if we can do this w/o evaluating the entire
   tree */


void multiplyBranchLengthsRecursively(TreeAln* traln, nodeptr p, double multiplier)
{  
  tree *tr = traln->getTr();
  double newZ = pow( traln->getBranchLength( p->number,0),multiplier); 
  // p->z[0] = p->back->z[0]
  traln->setBranchLengthSave(newZ, 0, p); 

  if(isTip(p->number, tr->mxtips))
    return; 

  nodeptr q = p->next; 
  while(p != q)
    {
      multiplyBranchLengthsRecursively(traln, q->back,multiplier); 
      q = q->next; 
    }
}


void applyTLMult(state *chain,  proposalFunction *pf )
{
  tree *tr = chain->traln->getTr(); 
  double multiplier = drawMultiplier(chain, pf->parameters.multiplier);
  pf->remembrance.multiplier = multiplier; 
  multiplyBranchLengthsRecursively(chain->traln , tr->start->back, multiplier); 
}


void resetTLMult(state *chain, proposalFunction *pf)
{
  tree *tr = chain->traln->getTr();
  multiplyBranchLengthsRecursively(chain->traln, tr->start->back, 1/pf->remembrance.multiplier); 
}
