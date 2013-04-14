#include "axml.h"
#include "bayes.h"
#include "adapters.h"
#include "randomness.h"
#include "output.h"
#include "path.h" 
#include "topology-utils.h"
#include "guidedMoves.h" 
#include "TreeAln.hpp"


/* 
   important TODO
   i am not sure, but treating branch lengths this way may seriously distort the lnls. 
   
   could be multiply the branches with a factor (2?), s.t. this gets better? 
   
 */



/**
    @brief finds the entry in the list 
    notice: inefficient, but who cares? 
 */ 
static insertList* getEntry(insertList *lnlList, branch b )
{
  for(insertList *iter = lnlList; iter ; iter = iter->next)
    {
      if(branchEqualUndirected(iter->b, b))
	return iter; 
    }
  return NULL; 
}


static void appendElem(insertList **list,  branch b , double lnl, double ratio, boolean isFirst  )
{
  insertList *elem = (insertList*)exa_calloc(1,sizeof(insertList)); 
  elem->next = *list; 
  *list = elem; 
  elem->lnl = lnl; 
  elem-> b = b; 
  elem->ratio = ratio;
  
  if(isFirst)
    elem->containedInFirst = TRUE; 
  else 
    elem->containedInSecond = TRUE; 
}

static void descendAndTestInsert(state *chain, branch pruneBranch, branch subtree, double ratio, insertList **lnlList, int depthLeft, boolean isFirst)
{ 
  tree *tr = chain->traln->getTr();

  /* insertList *la = *lnlList;  */
  if(depthLeft == 0 || isTip(pruneBranch.thisNode,tr->mxtips) )
    return; 
  
#if DEBUG_GUIDED_SPR > 1 
  printf("descendend to branch %d,%d\n", pruneBranch.thisNode, pruneBranch.thatNode);
#endif

  int
    numBranches = chain->traln->getNumBranches(); 
 
  nodeptr
    q = findNodeFromBranch(tr, pruneBranch), 
    p = findNodeFromBranch(tr, subtree); 

  /* first evaluation  */  
  branch descendBranch = constructBranch(q->next->back->number, q->number); 
  descendAndTestInsert(chain,descendBranch , subtree, 1  -  ratio, lnlList, depthLeft-1, isFirst ); 
  insertList *entry = NULL; 
  
  if( ( entry = getEntry(*lnlList, descendBranch) ) == NULL )
    {
      nodeptr 
	iP = q->next,    
	iP2 = q->next->back; 
      double zOld = iP->z[0],
	a,b;       
      divideBranchLengthsWithRatio(tr, zOld, ratio, &a, &b); 
      hookup(p->next,iP, &a, numBranches);   
      hookup(p->next->next,iP2, &b, numBranches); 

      newViewGenericWrapper(chain, p, FALSE); 
      evaluateGenericWrapper(chain, p, FALSE); 

      appendElem(lnlList, constructBranch(iP->number, iP2->number), tr->likelihood, ratio, isFirst); 
#if DEBUG_GUIDED_SPR > 1 
      printf("test insert: inserting %d into {%d,%d}(%g) => %.3f \n", p->number, iP->number, iP2->number, branchLengthToReal(tr, iP->z[0]), chain->tr->likelihood); 
#endif 
      /* remove the node again */
      hookup(iP, iP2, &zOld,numBranches); 
      newViewGenericWrapper(chain, iP, FALSE); /* OKAY  */
      newViewGenericWrapper(chain, iP2, FALSE);
      p->next->next->back = p->next->back = (nodeptr)NULL;       
    }
  else 
    {
      assert(NOT isFirst); 
      entry->containedInSecond = TRUE; 
      /* TODO assert this and that  */
    }


  /* second evaluation  */
  descendBranch = constructBranch(q->next->next->back->number, q->number); 
  descendAndTestInsert(chain, descendBranch, subtree, 1 - ratio, lnlList, depthLeft-1, isFirst ); 
  if( (entry = getEntry(*lnlList, descendBranch)) == NULL) /* we have not tested this branch  */
    {
      nodeptr iP = q->next->next, 
	iP2 = q->next->next->back; 

      double zOrig = iP->z[0],
	a,b; 
      divideBranchLengthsWithRatio(tr, zOrig, ratio, &a, &b); 

      hookup(p->next, iP, &a, numBranches); 
      hookup(p->next->next, iP2, &b, numBranches); 

      newViewGenericWrapper(chain,p,FALSE); 
      evaluateGenericWrapper(chain,p,FALSE); 

      appendElem(lnlList, constructBranch(iP->number, iP2->number), tr->likelihood, ratio, isFirst); 
#if DEBUG_GUIDED_SPR > 1 
      printf("test insert: inserting %d into {%d,%d}(%g) => %.3f\n", p->number, iP->number, iP2->number, branchLengthToReal(tr, iP->z[0]),  chain->tr->likelihood);   
#endif
      hookup(iP, iP2, &zOrig, numBranches); 
      newViewGenericWrapper(chain, iP, FALSE);
      newViewGenericWrapper(chain, iP2, FALSE);
      p->next->next->back = p->next->back = (nodeptr) NULL; 
    }
  else
    {
      assert(NOT isFirst); 
      entry->containedInSecond = TRUE; 
      /* TODO assert this and that  */
    }    
}


/**
   @brief transforms likelihoods into likelihood weights 
 */ 
void createWeights(tree  *tr, insertList **lnlList, boolean doFirst)
{
  insertList *la = *lnlList; 
  
  double best  = 1 ; 
  for(insertList *iter = la; iter; iter = iter->next)
    {
      if( ( doFirst && NOT iter->containedInFirst  ) 
	  || (NOT doFirst && NOT iter->containedInSecond) )
	continue; 

      if(best == 1 )
	best = iter->lnl; 
      else if(best < iter->lnl)
	best = iter->lnl; 
    }
  
  double sum = 0; 
  for(insertList *iter = la; iter; iter = iter->next)
    {
      if(doFirst && iter->containedInFirst)
	{
	  iter->weightInFirst = exp(iter->lnl - best ) + WEIGHT_EPS;
	  sum += iter->weightInFirst; 
	}
      else if(NOT doFirst && iter->containedInSecond)
	{
	  iter->weightInSecond = exp(iter->lnl - best ) + WEIGHT_EPS;
	  sum += iter->weightInSecond; 
	}
    }

  for(insertList *iter = la ; iter ; iter = iter->next)
    {
      if(doFirst && iter->containedInFirst) 
	iter->weightInFirst /= sum; 
      else if(NOT doFirst && iter->containedInSecond)
	iter->weightInSecond /= sum; 
    }  
}



static void testInsertWithRadius(state *chain, insertList **lnlList, branch subtree, branch pruneBranch, int radius, boolean isFirst)
{
  tree *tr = chain->traln->getTr();
  int numBranches = chain->traln->getNumBranches(); 
  assert(numBranches == 1 ); 

  double ratio = 0; 
  
  /* prune it  */
  {
    nodeptr 
      p = findNodeFromBranch(tr, subtree),
      q = p->next->back, 
      r = p->next->next->back; 

    ratio = getRatio(tr, q->z[0], r->z[0]);
    double tmp = combineBranchLengths(tr, q->z[0], r->z[0]);

    hookup(q,r, &tmp, numBranches); 
    p->next->back = p->next->next->back = (nodeptr)NULL;     
    newViewGenericWrapper(chain, q, FALSE); /* TODO no newview on r necessary?   */
  }

  descendAndTestInsert(chain, pruneBranch, subtree, ratio, lnlList,  radius, isFirst); 
  descendAndTestInsert(chain, invertBranch(pruneBranch), subtree, 1 - ratio, lnlList, radius, isFirst ); 

  /* reset the tree  */
  {
    nodeptr
      p = findNodeFromBranch(tr, subtree),
      q = findNodeFromBranch(tr, pruneBranch), 
      r = findNodeFromBranch(tr, invertBranch(pruneBranch)); 

    double z1, z2; 
    divideBranchLengthsWithRatio(tr, q->z[0], ratio, &z1, &z2); 
    hookup(p->next, q, &z1, numBranches); 
    hookup(p->next->next, r, &z2, numBranches);   
  }

  debug_checkTreeConsistency(chain->traln->getTr());
}



/**
   @brief samples a branch according to the likelihoods we determined
 */ 
insertList* sampleFromLnlArray(state *chain, insertList* lnlList)
{
  double
    r = drawRandDouble01(chain); 

  for(insertList *iter = lnlList; iter ; iter = iter->next)
    {
      double weight =  iter->weightInFirst; 
      if(r < weight)
	{
	  assert(iter->containedInFirst && NOT iter->containedInSecond); /* cannot be */
	  return iter; 
	}
      else 
	r -= weight; 
    }

  assert(0); 
  return NULL; 
}


static void freeList(insertList *lnlList)
{
  insertList *iter = lnlList; 
  while(iter) 
    {
      insertList *tmp = iter->next; 
      exa_free(iter); 
      iter = tmp; 
    }
}


void printInsertList(insertList *list)
{
  for(insertList *iter = list; iter ; iter = iter->next)
    {
      printf("{%d,%d}:lnl=%.2f, weights=[%g,%g]\n", iter->b.thisNode, iter->b.thatNode, 
	     iter->lnl, iter->weightInFirst, iter->weightInSecond); 
    }  
}



/**
   @brief a dummy evaluate for the guided radius spr move.
   
   actually this method should not do anything. Currently, we need
   it, to update the lnl arrays appropriatly. Pretty expensive.
*/ 
void evalGuidedSPR(state *chain, proposalFunction *pf)
{
  tree
    *tr = chain->traln->getTr(); 

  nodeptr p = findNodeFromBranch(tr, pf->remembrance.modifiedPath->content[1]);
  evaluateGenericWrapper(chain, p, FALSE );
}



void resetGuidedSPR(state *chain, proposalFunction *pf)
{
  tree *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches(); 
  assert(numBranches == 1); 

  path *rememPath = pf->remembrance.modifiedPath; 
  double ratio = pf->ratio ; 

  branch pruneBranch = rememPath->content[0],
    subtreeBranch = rememPath->content[1]; 

  nodeptr p = findNodeFromBranch(tr, subtreeBranch ); 
  double restoredZ = combineBranchLengths(tr, p->next->back->z[0], p->next->next->back->z[0]) ; 
  hookup(p->next->back, p->next->next->back, &restoredZ, numBranches); 
#if DEBUG_GUIDED_SPR > 1   
  printf("RESET: hooking up %d with %d (%g)\n", p->next->back->number, p->next->next->back->number, branchLengthToReal(tr, restoredZ)); 
#endif
  
  nodeptr q = findNodeFromBranch(tr, pruneBranch),
    r = findNodeFromBranch(tr, invertBranch(pruneBranch)); 
  double a,b; 
  divideBranchLengthsWithRatio(tr, q->z[0], ratio, &a,&b); 
  hookup(p->next, q, &a,numBranches); 
  hookup(p->next->next, r, &b,numBranches); 

#if DEBUG_GUIDED_SPR > 1   
  printf("RESET: hooking %d with %d,%d with %g,%g\tratio=%g\n", p->number, q->number, r->number, 
    branchLengthToReal(tr ,a),branchLengthToReal(tr, b), ratio); 
#endif
  
  debug_checkTreeConsistency(chain->traln->getTr());
#ifdef  DEBUG_SHOW_TREE 
  printf("RESET:\t");      
#endif
  debug_printTree(chain);
}



void debug_checkHooks(tree *tr)
{
#ifdef DEBUG_CHECK_TREE_CONSISTENCY
  for(int i = 1 ; i < 2 * tr->mxtips -2 ; ++i)
    {
      nodeptr p = tr->nodep[i]; 
      if(isTip(p->number, tr->mxtips))
	{
	  assert(p->back != NULL); 
	}
      else 
	{
	  assert(p->back != NULL 
		 && p->next->back != NULL
		 && p->next->next->back != NULL); 
	}
    }
#endif
}



void applyGuidedSPR(state *chain, proposalFunction *pf)
{  
#ifdef DEBUG_SHOW_TREE 
  printf("\n"); 
  double
    treeLength = getTreeLength(chain->tr, chain->tr->nodep[1]->back); 
  printf("tree length = %g\n", branchLengthToReal(chain->tr,treeLength)); 
  printf("BEFORE:\t"); 
#endif
  debug_printTree(chain);

  tree *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches(); 
  assert(numBranches == 1); 

  /* find subtree to prune  */
  boolean accepted = FALSE;  
  branch
    subtree, pruneBranch;   
  nodeptr p = NULL ; 
  double ratio = 0;
  
  while(NOT accepted) 
    {
      subtree = drawBranchUniform(chain);      
      p = findNodeFromBranch(tr, subtree); 
      accepted = NOT isTip(p->number, tr->mxtips)
	&& (NOT isTip(p->next->back->number,tr->mxtips ) && NOT isTip(p->next->next->back->number, tr->mxtips)); 
      if(accepted)
	{
	  pruneBranch.thisNode = p->next->back->number; 
	  pruneBranch.thatNode = p->next->next->back->number; 	  
	  ratio = getRatio(tr,  p->next->back->z[0] , p->next->next->back->z[0]);
	  /* printf("\nRATIO=%g\n", ratio) ; */
	}
    }


#if DEBUG_GUIDED_SPR > 0
  printf("picked subtree {%d,%d} to be pruned from {%d,%d} (%g,%g)\n", subtree.thisNode, subtree.thatNode, pruneBranch.thisNode, pruneBranch.thatNode, branchLengthToReal(tr, p->next->back->z[0]), branchLengthToReal(tr, p->next->next->back->z[0])); 
#endif


  insertList *lnlList = NULL; 
  testInsertWithRadius(chain, &lnlList,  subtree,  pruneBranch, pf->parameters.radius, TRUE ) ; 
  createWeights(tr,&lnlList, TRUE); 

  insertList *sampledElem = sampleFromLnlArray(chain, lnlList);
  branch sampledBranch = sampledElem->b; 

  /* prune subtree and insert into new position  */
  {
    nodeptr p = findNodeFromBranch(tr, subtree); 
    double z = combineBranchLengths(tr, p->next->back->z[0],   p->next->next->back->z[0]); 
    hookup(p->next->back, p->next->next->back, &z, numBranches); 
    
    nodeptr iP = findNodeFromBranch(tr, sampledBranch),
      iP2 = findNodeFromBranch(tr, invertBranch(sampledBranch)); 

    double a,b; 
    divideBranchLengthsWithRatio(tr, iP->z[0], sampledElem->ratio, &a,&b); 

#if DEBUG_GUIDED_SPR > 0 
    printf("inserting into %d into %d,%d(%g) with %g,%g\n", p->number, iP->number, iP2->number, branchLengthToReal(tr,iP->z[0]), branchLengthToReal(tr, a),branchLengthToReal(tr, b)); 
#endif

    hookup(p->next, iP, &a,numBranches); 
    hookup(p->next->next, iP2, &b, numBranches);     
    
  }
  
  testInsertWithRadius(chain, &lnlList, subtree, sampledBranch, pf->parameters.radius, FALSE); 
  createWeights(tr, &lnlList, FALSE); 



  /* do a uniform slide on the insertion branch: if we do not do that,
     we may get into a situation, where BLs become extremely short
     (when we only use guidedSPR+BLmult), because we have to accept
     the proposals, even with weird BLs.  */
  {
    nodeptr p = findNodeFromBranch(tr, subtree); 
    double sum = branchLengthToReal(tr, p->next->z[0] * p->next->next->z[0]);
    double onePart = drawRandDouble01(chain) * sum ; 
    double otherPart = sum - onePart; 
    onePart = branchLengthToInternal(tr, onePart); 
    otherPart = branchLengthToInternal(tr, otherPart); 
    hookup(p->next, p->next->back, &onePart, numBranches); 
    hookup(p->next->next, p->next->next->back, &otherPart, numBranches); 
    newViewGenericWrapper(chain, p, FALSE);
  }

#if DEBUG_GUIDED_SPR > 1  
  printInsertList(lnlList); 
#endif

  /* set the hastings */  
  insertList *entry = getEntry(lnlList, pruneBranch); 
  assert( sampledElem->containedInFirst && entry->containedInSecond  && chain->hastings == 1.); 
  assert(entry->weightInSecond != 0 &&  sampledElem->weightInFirst != 0); 
  chain->hastings = entry->weightInSecond / sampledElem->weightInFirst; 
#if DEBUG_GUIDED_SPR   > 1 
  printf("hastings: %g / %g => %g\n", entry->weightInSecond, sampledElem->weightInFirst , chain->hastings); 
#endif
  
  /* set this for the eval function */
  tr->likelihood = sampledElem->lnl;
  pf->ratio = ratio;  

  /* describe how to reset the move */
  path *rememPath = pf->remembrance.modifiedPath; 
  clearStack(rememPath); 
  pushStack(rememPath, pruneBranch); 
  pushStack(rememPath, subtree); 

  debug_checkHooks( tr); 
  freeList(lnlList);   
}
