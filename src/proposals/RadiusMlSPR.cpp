#include "RadiusMlSPR.hpp"
#include "eval.h"
// #include "output.h"


/* 
   important TODO
   i am not sure, but treating branch lengths this way may seriously distort the lnls. 
   
   could be multiply the branches with a factor (2?), s.t. this gets better? 
   
 */




RadiusMlSPR::RadiusMlSPR(  int _radius)
  : radius(_radius)
{
  this->name = "radiusMLspr" ;
  category = TOPOLOGY; 
  relativeWeight = 0.0;
  // ptype = GUIDED_SPR;   
}


static void divideBranchLengthsWithRatio(double orig, double ratio, double &resultA , double &resultB)
{
  resultA = pow(orig, (double)STRETCH_FACTOR * ratio) ; 
  resultB = pow(orig, (double)STRETCH_FACTOR * (1-ratio)); 

}



double getRatio(tree *tr, double a, double b )
{
  // sorry...not in the mood 
#if 1 
  return 1 ; 

#else 
  double realA = branchLengthToReal(tr, a ) , 
    realB = branchLengthToReal(tr,b); 
  
  double ratio = realA / (realA + realB ); 

  if(NOT (0 < ratio && ratio < 1.))
    {
      printf("\n\nWARNING: %g,%g (%g,%g) in ratio %g\n", a,b,branchLengthToReal(tr, a),branchLengthToReal(tr,b), ratio);
    }    
  return ratio ; 
#endif 
}


/**
    @brief finds the entry in the list 
    notice: inefficient, but who cares? 
 */ 
static insertList* getEntry(insertList *lnlList, Branch b )
{
  for(insertList *iter = lnlList; iter ; iter = iter->next)
    {
      if(iter->b.equalsUndirected( b))
	return iter; 
    }
  return NULL; 
}


static void appendElem(insertList **list,  Branch b , double lnl, double ratio, boolean isFirst  )
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

static void descendAndTestInsert(TreeAln& traln, Branch pruneBranch, Branch subtree, double ratio, insertList **lnlList, int depthLeft, boolean isFirst)
{ 
  tree *tr = traln.getTr();

  if( depthLeft == 0 || isTip(pruneBranch.getPrimNode(),tr->mxtips) )
    return; 
  
#if DEBUG_GUIDED_SPR > 1 
  printf("descendend to branch %d,%d\n", pruneBranch.thisNode, pruneBranch.thatNode);
#endif

  nodeptr
    q = pruneBranch.findNodePtr(traln ), 
    p = subtree.findNodePtr(traln ); 

  /* first evaluation  */  
  Branch descendBranch(q->next->back->number, q->number); 
  descendAndTestInsert(traln,descendBranch , subtree, 1  -  ratio, lnlList, depthLeft-1, isFirst ); 
  insertList *entry = NULL; 
  
  if( ( entry = getEntry(*lnlList, descendBranch) ) == NULL )
    {
      nodeptr 
	iP = q->next,    
	iP2 = q->next->back; 
      double zOld = traln.getBranchLength(iP, 0),
	a,b; 
      divideBranchLengthsWithRatio( zOld, ratio, a, b); 
      traln.clipNode(p->next,iP, a);   
      traln.clipNode(p->next->next,iP2, b); 
      traln.setBranchLengthBounded(a, 0, p->next );
      traln.setBranchLengthBounded(b, 0, p->next->next );

      newViewGenericWrapper(traln, p, FALSE); 
#if 0 
      evaluateGenericWrapper(traln, p, FALSE ); 
#endif

      appendElem(lnlList, Branch(iP->number, iP2->number), tr->likelihood, ratio, isFirst); 
#if DEBUG_GUIDED_SPR > 1 
      printf("test insert: inserting %d into {%d,%d}(%g) => %.3f \n", p->number, iP->number, iP2->number, branchLengthToReal(tr, traln.getBranchLength(  iP,0)), tr->likelihood); 
#endif 
      /* remove the node again */
      traln.clipNode(iP, iP2, zOld); 
      newViewGenericWrapper(traln, iP, FALSE); /* OKAY  */
      newViewGenericWrapper(traln, iP2, FALSE);
      p->next->next->back = p->next->back = (nodeptr)NULL;       
    }
  else 
    {
      assert(NOT isFirst); 
      entry->containedInSecond = TRUE; 
      /* TODO assert this and that  */
    }


  /* second evaluation  */
  descendBranch = Branch(q->next->next->back->number, q->number); 
  descendAndTestInsert(traln, descendBranch, subtree, 1 - ratio, lnlList, depthLeft-1, isFirst ); 
  if( (entry = getEntry(*lnlList, descendBranch)) == NULL) /* we have not tested this branch  */
    {
      nodeptr iP = q->next->next, 
	iP2 = q->next->next->back; 

      double zOrig = traln.getBranchLength(iP,0),
	a,b; 
      divideBranchLengthsWithRatio( zOrig, ratio, a, b); 

      traln.clipNode(p->next, iP, a); 
      traln.clipNode(p->next->next, iP2, b); 
      traln.setBranchLengthBounded(a,0,p->next);
      traln.setBranchLengthBounded(b,0,p->next->next);

      newViewGenericWrapper(traln,p,FALSE); 
#if 0 
      evaluateGenericWrapper(traln,p,FALSE  ); 
#endif

      appendElem(lnlList, Branch(iP->number, iP2->number), tr->likelihood, ratio, isFirst); 
#if DEBUG_GUIDED_SPR > 1 
      printf("test insert: inserting %d into {%d,%d}(%g) => %.3f\n", p->number, iP->number, iP2->number, branchLengthToReal(tr, traln.getBranchLength(iP->number, 0)),  tr->likelihood);   
#endif
      traln.clipNode(iP, iP2, zOrig); 
      newViewGenericWrapper(traln, iP, FALSE);
      newViewGenericWrapper(traln, iP2, FALSE);
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
static  void createWeights(tree  *tr, insertList **lnlList, boolean doFirst)
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




static void testInsertWithRadius(TreeAln &traln, insertList **lnlList, Branch subtree, Branch pruneBranch, int radius, boolean isFirst)
{
#if 0 
  tree *tr = traln.getTr();
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 

  double ratio = 0; 
  
  /* prune it  */
  {
    nodeptr 
      p = subtree.findNodePtr(traln ),
      q = p->next->back, 
      r = p->next->next->back; 

    ratio = getRatio(tr, traln.getBranchLength(q,0), traln.getBranchLength(r,0));
    double tmp = combineBranchLengths(tr, traln.getBranchLength(q,0), traln.getBranchLength(r,0));

    traln.clipNode(q,r, tmp); 
    p->next->back = p->next->next->back = (nodeptr)NULL;     
    newViewGenericWrapper(traln, q, FALSE); /* TODO no newview on r necessary?   */
  }

  descendAndTestInsert(traln, pruneBranch, subtree, ratio, lnlList,  radius, isFirst); 
  descendAndTestInsert(traln, invertBranch(pruneBranch), subtree, 1 - ratio, lnlList, radius, isFirst ); 

  /* reset the tree  */
  {
    nodeptr
      p = findNodeFromBranch(tr, subtree),
      q = findNodeFromBranch(tr, pruneBranch), 
      r = findNodeFromBranch(tr, invertBranch(pruneBranch)); 

    double z1, z2; 
    divideBranchLengthsWithRatio(tr, traln.getBranchLength(q,0), ratio, &z1, &z2); 
    traln.clipNode(p->next, q, z1); 
    traln.clipNode(p->next->next, r, z2);   
    traln.setBranchLengthBounded(z1,0,p->next);
    traln.setBranchLengthBounded(z2,0,p->next->next); 
  }

  debug_checkTreeConsistency(traln);
#endif
}



/**
   @brief samples a branch according to the likelihoods we determined
 */ 
static insertList* sampleFromLnlArray(TreeAln& traln, Randomness &rand, insertList* lnlList)
{
  double
    r = rand.drawRandDouble01(); 

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


static void debug_checkHooks(tree *tr)
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



void RadiusMlSPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
#if 0 
  debug_printTree(traln);

  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1); 

  /* find subtree to prune  */
  boolean accepted = FALSE;  
  Branch
    subtree, pruneBranch;   
  nodeptr p = NULL ; 
  double ratio = 0;
  
  while(NOT accepted) 
    {
      subtree = traln.drawBranchUniform(rand); 
      p = subtree.findNodePtr(traln ); 
      accepted = NOT isTip(p->number, tr->mxtips)
	&& (NOT isTip(p->next->back->number,tr->mxtips ) && NOT isTip(p->next->next->back->number, tr->mxtips)); 
      if(accepted)
	{
	  pruneBranch.setPrimNode( p->next->back->number); 
	  pruneBranch.setSecNode(p->next->next->back->number); 
	  ratio = getRatio(tr,  traln.getBranchLength(p->next->back, 0 ) , traln.getBranchLength( p->next->next->back,0));
	  /* printf("\nRATIO=%g\n", ratio) ; */
	}
    }


#if DEBUG_GUIDED_SPR > 0
  printf("picked subtree {%d,%d} to be pruned from {%d,%d} (%g,%g)\n", subtree.thisNode, subtree.thatNode, pruneBranch.thisNode, pruneBranch.thatNode, branchLengthToReal(tr, traln.getBranchLength( p->next->back->number,0)), branchLengthToReal(tr, traln.getBranchLength( p->next->next->back->number,0))); 
#endif


  insertList *lnlList = NULL; 
  testInsertWithRadius(traln, &lnlList,  subtree,  pruneBranch, radius, TRUE ) ; 
  createWeights(tr,&lnlList, TRUE); 

  insertList *sampledElem = sampleFromLnlArray(traln, rand, lnlList);
  Branch sampledBranch = sampledElem->b; 

  /* prune subtree and insert into new position  */
  {
    nodeptr p = subtree.findNodePtr(traln ); 
    double z = combineBranchLengths(tr, traln.getBranchLength( p->next->back,0),   traln.getBranchLength( p->next->next->back,0)); 
    traln.clipNode(p->next->back, p->next->next->back, z); 
    
    nodeptr iP = findNodeFromBranch(tr, sampledBranch),
      iP2 = findNodeFromBranch(tr, invertBranch(sampledBranch)); 

    double a,b; 
    divideBranchLengthsWithRatio(tr, traln.getBranchLength( iP,0), sampledElem->ratio, &a,&b); 

#if DEBUG_GUIDED_SPR > 0 
    printf("inserting into %d into %d,%d(%g) with %g,%g\n", p->number, iP->number, iP2->number, branchLengthToReal(tr, traln.getBranchLength( iP->number,0)), branchLengthToReal(tr, a),branchLengthToReal(tr, b)); 
#endif

    traln.clipNode(p->next, iP, a); 
    traln.setBranchLengthBounded(a,0,p->next);
    traln.clipNode(p->next->next, iP2, b);     
    traln.setBranchLengthBounded(b,0,p->next->next);    
  }
  
  testInsertWithRadius(traln, &lnlList, subtree, sampledBranch, radius, FALSE); 
  createWeights(tr, &lnlList, FALSE); 

  // TODO prior 
  assert(0); 

  /* do a uniform slide on the insertion branch: if we do not do that,
     we may get into a situation, where BLs become extremely short
     (when we only use guidedSPR+BLmult), because we have to accept
     the proposals, even with weird BLs.  */
  {
    nodeptr p = findNodeFromBranch(tr, subtree); 
    double sum = branchLengthToReal(tr, traln.getBranchLength( p->next,0) * traln.getBranchLength( p->next->next,0));
    double onePart = rand.drawRandDouble01() * sum ; 
    double otherPart = sum - onePart; 
    onePart = branchLengthToInternal(tr, onePart); 
    otherPart = branchLengthToInternal(tr, otherPart); 
    traln.clipNode(p->next, p->next->back, onePart); 
    traln.clipNode(p->next->next, p->next->next->back, otherPart); 
    newViewGenericWrapper(traln, p, FALSE);
  }

#if DEBUG_GUIDED_SPR > 1  
  printInsertList(lnlList); 
#endif

  /* set the hastings */  
  insertList *entry = getEntry(lnlList, pruneBranch); 
  assert( sampledElem->containedInFirst && entry->containedInSecond); 
  assert(entry->weightInSecond != 0 &&  sampledElem->weightInFirst != 0); 
  
  updateHastings(hastings,entry->weightInSecond / sampledElem->weightInFirst , name ); 
#if DEBUG_GUIDED_SPR   > 1 
  printf("ln-hastings: %g / %g => %g\n", entry->weightInSecond, sampledElem->weightInFirst , chain->hastings); 
#endif
  
  /* set this for the eval function */
  tr->likelihood = sampledElem->lnl;
  this->ratio = ratio;  

  /* describe how to reset the move */
  // Path *rememPath = pf->remembrance.modifiedPath; 
  path.clear(); 
  path.append( pruneBranch); 
  path.append( subtree); 

  debug_checkHooks( tr); 
  freeList(lnlList);   
#endif
} 


/**
   @brief a dummy evaluate for the guided radius spr move.
   
   actually this method should not do anything. Currently, we need
   it, to update the lnl arrays appropriatly. Pretty expensive.
*/ 
void RadiusMlSPR::evaluateProposal(  LikelihoodEvaluatorPtr &evaluator, TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = path.at(1).findNodePtr(traln );
#if 0 
  evaluateGenericWrapper(traln, p, FALSE );
#endif
} 


void RadiusMlSPR::resetState(TreeAln &traln, PriorBelief &prior) 
{
#if 0 
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1); 
  
  Branch pruneBranch = path.at(0),
    subtreeBranch = path.at(1); 

  nodeptr p = findNodeFromBranch(tr, subtreeBranch ); 
  double restoredZ = combineBranchLengths(tr, traln.getBranchLength(p->next->back, 0), traln.getBranchLength(p->next->next->back, 0)) ; 
  traln.clipNode(p->next->back, p->next->next->back, restoredZ); 
#if DEBUG_GUIDED_SPR > 1   
  printf("RESET: hooking up %d with %d (%g)\n", p->next->back->number, p->next->next->back->number, branchLengthToReal(tr, restoredZ)); 
#endif
  
  nodeptr q = findNodeFromBranch(tr, pruneBranch),
    r = findNodeFromBranch(tr, invertBranch(pruneBranch)); 
  double a,b; 
  divideBranchLengthsWithRatio(tr, traln.getBranchLength(q, 0 ), ratio, &a,&b); 
  traln.clipNode(p->next, q, a); 
  traln.clipNode(p->next->next, r, b); 

#if DEBUG_GUIDED_SPR > 1   
  printf("RESET: hooking %d with %d,%d with %g,%g\tratio=%g\n", p->number, q->number, r->number, 
    branchLengthToReal(tr ,a),branchLengthToReal(tr, b), ratio); 
#endif
  
  debug_checkTreeConsistency(traln);
#ifdef  DEBUG_SHOW_TREE 
  printf("RESET:\t");      
#endif
  debug_printTree(traln);
#endif
} 



AbstractProposal* RadiusMlSPR::clone() const
{
  // return new RadiusMlSPR( radius);
  return new RadiusMlSPR( *this );
}

