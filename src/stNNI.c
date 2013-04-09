

#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "path.h"
#include "output.h"

void apply_st_nni(state *chain, proposalFunction *pf)
{  
  tree *tr = chain->tr; 
  int numBranches = getNumBranches(tr);
  assert(numBranches == 1); 
  branch b = drawInnerBranchUniform(chain); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  branch toBePruned = constructBranch(p->next->back->number, p->next->number ); 
  
  branch switchingBranch; 
  double r = drawRandDouble01(chain);
  if(r < 0.5)
    switchingBranch = constructBranch(p->back->next->back->number, p->back->next->number); 
  else 
    switchingBranch = constructBranch(p->back->next->next->back->number, p->back->next->next->number); 
  
  path *rememStack = pf->remembrance.modifiedPath; 
  clearStack(rememStack); 
  pushStack(rememStack, toBePruned);
  pushStack(rememStack, switchingBranch);

  /* switch the branches  */
  {
    nodeptr q = findNodeFromBranch(tr, toBePruned),
      qBack = q->back; 
    
    nodeptr r = findNodeFromBranch(tr, switchingBranch),
      rBack = r->back; 

    hookup(r, qBack, r->z, numBranches);
    hookup(q, rBack, q->z, numBranches);    
  }

  debug_checkTreeConsistency(chain);
  /* TODO maybe multiply as well */
  
}

void eval_st_nni(state *chain, proposalFunction *pf )
{
  tree *tr = chain->tr; 

  branch b1 = pf->remembrance.modifiedPath->content[0],
    b2  = pf->remembrance.modifiedPath->content[1]; 
  
  branch exchangeBranch = constructBranch(b1.thatNode, b2.thatNode);
  
  nodeptr p = findNodeFromBranch(tr, exchangeBranch),
    q = findNodeFromBranch(tr, invertBranch(exchangeBranch)); 

  exa_newViewGeneric(chain, p, FALSE); 
  exa_newViewGeneric(chain, q, FALSE); 

  evaluateGenericWrapper(chain,p,FALSE); 
  /* printf("lnl = %g\n", tr->likelihood);  */
}

void reset_st_nni(state *chain, proposalFunction *pf)
{
  tree
    *tr = chain->tr; 
  int numBranches = getNumBranches(tr); 
  
  path *rStack = pf->remembrance.modifiedPath; 
  branch a = popStack(rStack),	/* switchingBranch */
    b = popStack(rStack); 	/* toBePruned */

  swpInt(&a.thatNode, &b.thatNode); 

  /* reset the branches */
  {
    nodeptr q = findNodeFromBranch(tr, a),
      qBack = q->back; 
    nodeptr r = findNodeFromBranch(tr, b), 
      rBack = r->back; 

    hookup(r, qBack, r->z, numBranches); 
    hookup(q, rBack, q->z, numBranches) ; 
  }

  debug_checkTreeConsistency(chain);

}
