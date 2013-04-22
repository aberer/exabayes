
#include "path.h"

#include "TreeAln.hpp"



/**
   @brief is node an outer node in the path?   
 */
boolean isOuterNode(int node, path *aPath)
{
  boolean occured = FALSE; 
  for(int i = 0;  i < aPath->index; ++i)
    {
      if(nodeIsInBranch(node, aPath->content[i]))
	{
	  if(occured)
	    return TRUE; 
	  else 
	    occured = TRUE; 
	}	
    }

  return FALSE; 
}


/* TODO should be static  */
void pushToStackIfNovel(path *s, branch b, int numTip)
{  
  assert(s->index >= 2); 

  branch bPrev = peekStack(s); 
  if(s->index == 2)
    {
      if(NOT branchEqualUndirected(b, bPrev))
	pushStack(s,b); 
    }
  else 
    {      
      if(branchEqualUndirected(b,bPrev))
	{
	  popStack(s); 
	}
      else if(isTipBranch(bPrev, numTip))
	{
	  popStack(s); 
	  pushStack(s,b); 
	}
      else 
	pushStack(s, b ); 
    }
}



boolean nodeIsOnPath(int node, path *aPath)
{
  for(int i = 0; i < aPath->index; ++i)
    {
      if(nodeIsInBranch(node, aPath->content[i]))
	return TRUE;       
    }

  return FALSE; 
}


/**
    @brief saves all branch lengths along the path in s. 
 */ 
void saveBranchLengthsPath(state *chain, path *s)
{
  TreeAln *traln = chain->traln; 
  tree *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches() ;
  
  for(int i = 0; i < s->index; ++i)
    {
      branch *bPtr = s->content + i ; 
      nodeptr p = findNodeFromBranch(tr, *bPtr); 

      for(int j = 0; j < numBranches; ++j)
	bPtr->length[j] = traln->getBranchLength( p,0); 
    }
}


/* draws a random path => would be practical */
void generateRandomPath(state *chain ,path *s, double stopProp)
{
  assert(0); 
  
  
}

static void debug_assertPathExists(tree *tr, path *s)
{
#ifdef DEBUG_CHECK_TREE_CONSISTENCY
  for(int i = 0; i < s->index; ++i)
    assert(branchExists(tr, s->content[i])); 
#endif
}



static void multiplyBranch(state *chain, branch b, double parameter, double *hastings)
{
  TreeAln *traln = chain->traln; 
  tree *tr = chain->traln->getTr(); 
  int numBranches = chain->traln->getNumBranches();
  nodeptr  p = findNodeFromBranch(tr, b); 
  double multiplier = drawMultiplier(chain, parameter); 
  double newZ = branchLengthToInternal(tr, multiplier * branchLengthToReal(tr, traln->getBranchLength( p,0))); 

  *hastings *= multiplier; 
  hookup(p,p->back, &newZ, numBranches);   
}


/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void multiplyAlongBranchESPR(state *chain, path *s, double multi )
{
  /* printf("multiplyAlongBranchESPR \n");  */
  tree *tr = chain->traln->getTr(); 
  assert(s->index >= 2); 
  int numBranches = chain->traln->getNumBranches(); 
  assert(numBranches == 1 ); 

  /* first two branches: consider that subtree has been pruned */
  branch b = getThirdBranch(tr, s->content[0], s->content[1]); 
  int sTNode = b.thisNode; 

  branch firstBranch = constructBranch( getOtherNode(sTNode, s->content[0]) ,
					getOtherNode(sTNode, s->content[1]));   
  double hastings = 1; 

  /* TODO hack: using the branch length multiplier parameter here  */
  /* proposalFunction *pf = NULL;  */
  /* for(int i = 0; i < chain->numProposals; ++i) */
  /*   { */
  /*     proposalFunction *tmp = chain->proposals[i];  */
  /*     if(tmp->ptype == BRANCH_LENGTHS_MULTIPLIER) */
  /* 	pf= tmp ;  */
  /*   } */
  /* if(NOT pf) */
  /*   { */
  /*     /\* printf("sorry, please activate the BL multiplier proposal.\n");  *\/ */
  /*     /\* assert(pf);  *\/ */
  /*   } */
  /* double parameter = pf->param2.multiplier ;  */
  
  
  multiplyBranch(chain, firstBranch, multi, &hastings); 
  
  /* treat all branches except the first 2 and the last one */
  for(int pathPos = 2; pathPos < s->index-1; ++pathPos)
    {
      branch b = s->content[pathPos]; 
      multiplyBranch(chain, b, multi, &hastings); 
    }

  /* treat last two paths (only one in representation) */
  b = constructBranch(sTNode, s->content[s->index-1].thisNode); 
  multiplyBranch(chain, b,multi, &hastings); 

  b = constructBranch(sTNode, s->content[s->index-1].thatNode); 
  multiplyBranch(chain,b,multi, &hastings);   
  chain->hastings *= hastings; 
}



/**
   @brief draws a random set of branches in the tree that constitute a
   path
   
   This function employs the eSPR strategy 
   
   @param s -- the result : the first two branches in the stack define
   the root of the pruned subtree
 */
void drawPathForESPR(state *chain,path *s, double stopProp )
{
  TreeAln *traln = chain->traln; 

  assert(0 < stopProp && stopProp < 1.0 ); 

  assert(stackIsEmpty(s)); 
  tree *tr  = chain->traln->getTr();   

  branch start; 
  nodeptr p,q,r; 
  do 
    {
      start  = drawSubtreeUniform(chain);
  
      p = findNodeFromBranch(tr, start);
      q = p->next->back; 
      r = p->next->next->back;
    } while(isTip(q->number,tr->mxtips) && isTip(r->number, tr->mxtips)); 

  /* TODO assert that the pruned subtree is not too big (s.t. we cannot find any insertion points)  */


  /* save branches and prune */
  double zqr[NUM_BRANCHES]; 
  for(int i = 0; i < chain->traln->getNumBranches(); ++i)
    zqr[i] = traln->getBranchLength( r,0); 
  hookup(q,r, q->z, chain->traln->getNumBranches());
  p->next->back = p->next->next->back = (nodeptr)NULL; 

  pushStack(s, start); 
  pushStack(s,constructBranch(q->number, r->number)); 

  nodeptr currentNode = drawRandDouble01(chain) ? q : r; 
  boolean accepted = FALSE;   
  /* printf("%d\n", currentNode->number); */
  while(NOT accepted)
    {
      nodeptr n = 
	drawRandDouble01(chain)  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      pushToStackIfNovel(s, constructBranch(currentNode->number, n->number),tr->mxtips); 
      
      /* printf("pushing %d,%d\n", currentNode->number, n->number);  */

      currentNode = n; 
      
      accepted = drawRandDouble01(chain) < stopProp && stackLength(s) > 2 ; 	
      /* printStack(s); */
    }

  /* undo changes to the tree  */
  hookup(p->next,q, q->z, chain->traln->getNumBranches()); 
  hookup(p->next->next,r,zqr, chain->traln->getNumBranches());   

  /* now correct  */
  if(nodeIsInBranch(s->content[1].thisNode, s->content[2] ))    
    s->content[1].thatNode = p->number; 
  else if(nodeIsInBranch(s->content[1].thatNode, s->content[2] ))
    s->content[1].thisNode = p->number; 
  else 
    assert(0); 

  /* correct the incorrectly set first branch in the path */
  s->content[0] = getThirdBranch(tr, s->content[0], s->content[1]); 
  /* printStack(s); */

  debug_assertPathExists(tr, s); 
}
