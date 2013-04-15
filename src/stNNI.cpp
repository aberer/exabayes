#include "axml.h"
#include "bayes.h"
#include "randomness.h"
#include "path.h"
#include "output.h"
#include "misc-utils.h"
#include "TreeAln.hpp"

// #define DEBUG_INFO 


void apply_st_nni(state *chain, proposalFunction *pf)
{  
  tree *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches();
  assert(numBranches == 1); 
  branch b = drawInnerBranchUniform(chain); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  branch toBePruned = constructBranch(p->next->back->number, p->number ); 
  
  branch switchingBranch; 
  double r = drawRandDouble01(chain);
  double bl = 0; 
  if(r < 0.5)
    {
      switchingBranch = constructBranch(p->back->next->back->number, p->back->number); 
      bl = p->back->next->back->z[0]; 
    }
  else 
    {
      switchingBranch = constructBranch(p->back->next->next->back->number, p->back->number); 
      bl = p->back->next->next->back->z[0]; 
    }

#ifdef DEBUG_INFO
  cout << "stNNI: switching " <<  toBePruned <<  " with "  << switchingBranch << "\t common branch " << b << endl; 
#endif

  toBePruned.length[0] = p->next->back->z[0]; 
  b.length[0] = p->z[0]; 
  switchingBranch.length[0] = bl; 
  
  path *rememStack = pf->remembrance.modifiedPath; 
  clearStack(rememStack); 
  pushStack(rememStack, toBePruned);
  pushStack(rememStack, b); 
  pushStack(rememStack, switchingBranch);

  /* switch the branches  */
  {
    nodeptr q = findNodeFromBranch(tr, toBePruned),
      qBack = q->back; 
    
    nodeptr r = findNodeFromBranch(tr, switchingBranch),
      rBack = r->back; 
    
    
#ifdef NNI_MULTIPLY_BL

    /* DIRTY */
    double m1 = drawMultiplier(chain, pf->parameters.multiplier),
      m2 =  drawMultiplier(chain, pf->parameters.multiplier),
      m3 =  drawMultiplier(chain, pf->parameters.multiplier);
      
    p->z[0] = pow(p->z[0],m1); 
    r->z[0] = pow(r->z[0],m2); 
    q->z[0] = pow(q->z[0],m3);     
#endif

    hookup(p,p->back, p->z, numBranches); 
    hookup(r, qBack, r->z, numBranches);
    hookup(q, rBack, q->z, numBranches);    

  }
  


  debug_checkTreeConsistency(chain->traln->getTr());
  /* TODO maybe multiply as well */
  
}

void eval_st_nni(state *chain, proposalFunction *pf )
{
  tree *tr = chain->traln->getTr(); 

  branch b1 = pf->remembrance.modifiedPath->content[0], // tobepruned 
    b2  = pf->remembrance.modifiedPath->content[2];	// switchingbranch 
  
  branch exchangeBranch = constructBranch(b1.thatNode, b2.thatNode);
  
  nodeptr p = findNodeFromBranch(tr, exchangeBranch),
    q = findNodeFromBranch(tr, invertBranch(exchangeBranch)); 

  newViewGenericWrapper(chain, p, FALSE); 
  newViewGenericWrapper(chain, q, FALSE); 

  evaluateGenericWrapper(chain,p,FALSE ); 
}



void reset_st_nni(state *chain, proposalFunction *pf)
{
  tree
    *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches(); 
  
  path *rStack = pf->remembrance.modifiedPath; 
  branch a = popStack(rStack),	/* switchingBranch */
    chosenBranch = popStack(rStack), /* inner branch   */
    b = popStack(rStack); 	/* toBePruned */

  swpInt(&a.thatNode, &b.thatNode); 

  /* reset the branches */
  {
    nodeptr q = findNodeFromBranch(tr, a),
      qBack = q->back; 
    nodeptr r = findNodeFromBranch(tr, b), 
      rBack = r->back; 
    
    nodeptr between = findNodeFromBranch(tr, chosenBranch ); 
    
    hookup(between, between->back, chosenBranch.length, numBranches); 
    hookup(r, qBack, b.length, numBranches);  /* r->z */
    hookup(q, rBack, a.length, numBranches) ; /* q->z */
  }

  debug_checkTreeConsistency(chain->traln->getTr());
}
