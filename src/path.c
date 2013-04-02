
#include "path.h"



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
  tree *tr = chain->tr; 
  int numBranches = getNumBranches(tr) ;
  
  for(int i = 0; i < s->index; ++i)
    {
      branch *bPtr = s->content + i ; 
      nodeptr p = findNodeFromBranch(tr, *bPtr); 

      for(int j = 0; j < numBranches; ++j)
	bPtr->length[j] = p->z[j]; 
    }
}


/* draws a random path => would be practical */
void generateRandomPath(state *chain ,path *s, double stopProp)
{
  assert(0); 
  
  
}

static void debug_assertPathExists(tree *tr, path *s)
{
#ifdef DEBUG_CHECK_TOPO_CONSISTENCY
  for(int i = 0; i < s->index; ++i)
    assert(branchExists(tr, s->content[i])); 
#endif
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
  assert(0 < stopProp && stopProp < 1.0 ); 

  assert(stackIsEmpty(s)); 
  tree *tr  = chain->tr;   

  branch start; 
  nodeptr p,q,r; 
  do 
    {
      start  = drawSubtreeUniform(chain);
  
      p = findNodeFromBranch(tr, start);
      q = p->next->back; 
      r = p->next->next->back;
      /* printf("\nchose subtree for pruning {%d,%d} => hooked to %d (%.2f),%d (%.2f)\n", start.thisNode, start.thatNode, p->next->back->number, p->next->back->z[0], p->next->next->back->number, p->next->next->back->z[0]); */
    } while(isTip(q->number,tr->mxtips) && isTip(r->number, tr->mxtips)); 

  /* TODO assert that the pruned subtree is not too big (s.t. we cannot find any insertion points)  */


  /* save branches and prune */
  double zqr[NUM_BRANCHES]; 
  for(int i = 0; i < getNumBranches(tr); ++i)
    zqr[i] = r->z[i]; 
  hookup(q,r, q->z, getNumBranches(tr));
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
  hookup(p->next,q,q->z, getNumBranches(tr)); 
  hookup(p->next->next,r,zqr, getNumBranches(tr));   

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




