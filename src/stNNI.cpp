#include "axml.h"
#include "bayes.h"
#include "Path.hpp"
#include "output.h"
#include "misc-utils.h"
#include "TreeAln.hpp"
#include "Randomness.hpp"

// #define DEBUG_INFO 


void apply_st_nni(Chain *chain, proposalFunction *pf)
{  
  TreeAln *traln = chain->traln; 
  tree *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches();
  assert(numBranches == 1); 
  branch b = chain->getChainRand()->drawInnerBranchUniform(*traln); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  branch toBePruned = constructBranch(p->next->back->number, p->number ); 
  
  branch switchingBranch; 
  double r = chain->getChainRand()->drawRandDouble01();
  double bl = 0; 
  if(r < 0.5)
    {
      switchingBranch = constructBranch(p->back->next->back->number, p->back->number); 
      bl = traln->getBranchLength( p->back->next->back,0); 
    }
  else 
    {
      switchingBranch = constructBranch(p->back->next->next->back->number, p->back->number); 
      bl = traln->getBranchLength( p->back->next->next->back,0); 
    }

#ifdef DEBUG_INFO
  cout << "stNNI: switching " <<  toBePruned <<  " with "  << switchingBranch << "\t common branch " << b << endl; 
#endif

  toBePruned.length[0] = traln->getBranchLength( p->next->back,0); 
  b.length[0] = traln->getBranchLength( p,0); 
  switchingBranch.length[0] = bl; 
  
  Path *rememStack = pf->remembrance.modifiedPath; 
  rememStack->clearStack(); 
  rememStack->pushStack( toBePruned);
  rememStack->pushStack( b); 
  rememStack->pushStack( switchingBranch);

  /* switch the branches  */
  {
    nodeptr q = findNodeFromBranch(tr, toBePruned),
      qBack = q->back; 
    
    nodeptr r = findNodeFromBranch(tr, switchingBranch),
      rBack = r->back; 
    
    
#ifdef NNI_MULTIPLY_BL

    /* DIRTY */
    double m1 = chain->getChainRand()->drawMultiplier( pf->parameters.multiplier),
      m2 =  chain->getChainRand()->drawMultiplier(pf->parameters.multiplier),
      m3 =  chain->getChainRand()->drawMultiplier(pf->parameters.multiplier);
      

    double old1 = traln->getBranchLength( p,0),
      old2 = traln->getBranchLength( r,0), 
      old3 = traln->getBranchLength( q,0),
      new1 = pow( old1,m1),
      new2 = pow( old2 ,m2),
      new3 = pow( old3 ,m3); 

#ifdef PRINT_MULT
    printf("%f*%f=%f\t%f*%f=%f\t%f*%f=%f\n", old1, m1, new1, old2, m2, new2, old3,m3,new3) ;
#endif

    traln->setBranchLengthSave( new1 ,0,p); 
    traln->setBranchLengthSave( new2 ,0,r); 
    traln->setBranchLengthSave( new3 ,0,q); 
    
#endif

    hookup(p,p->back, p->z, numBranches); 
    hookup(r, qBack, r->z, numBranches);
    hookup(q, rBack, q->z, numBranches);    

  }
  


  debug_checkTreeConsistency(chain->traln->getTr());
  /* TODO maybe multiply as well */
  
}



void eval_st_nni(Chain *chain, proposalFunction *pf )
{
  tree *tr = chain->traln->getTr(); 

  branch b1 = pf->remembrance.modifiedPath->at(0), // tobepruned 
    b2  = pf->remembrance.modifiedPath->at(2);	// switchingbranch 
  
  branch exchangeBranch = constructBranch(b1.thatNode, b2.thatNode);
  
  nodeptr p = findNodeFromBranch(tr, exchangeBranch),
    q = findNodeFromBranch(tr, invertBranch(exchangeBranch)); 

  newViewGenericWrapper(chain, p, FALSE); 
  newViewGenericWrapper(chain, q, FALSE); 

  evaluateGenericWrapper(chain,p,FALSE ); 
}



void reset_st_nni(Chain *chain, proposalFunction *pf)
{
  tree
    *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches(); 
  
  Path *rStack = pf->remembrance.modifiedPath; 
  branch a = rStack->popStack(),	/* switchingBranch */
    chosenBranch = rStack->popStack(), /* inner branch   */
    b = rStack->popStack(); 	/* toBePruned */

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

