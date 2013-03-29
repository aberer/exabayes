/**
   @file branch.c
   
   @brief Everything that deals with branches.  
 */

#include "axml.h"
#include "bayes.h"
#include "branch.h" 
#include "adapters.h"
#include "output.h"  


static nodeptr findEmptyNodePtr(nodeptr ptr)
{
  if(ptr->back == NULL)
    return ptr;
  else if (ptr->next->back == NULL)
    return ptr->next;
  else if(ptr->next->next->back == NULL)
    return ptr->next->next;
  else
    {
      assert(0);
      return NULL;
    }
}

/**
   @brief Inserts a node into a branch for given node ids. 

   This could replace with the other insert function in this file, if
   we decide to stick to the node-id scheme.

   Notice that, two next pointers must be NULL, for the
   insertNode. This is the case by default, if you called
   pruneNodeFromNodes before.

   @param insertBranch -- the branch we insert (thisNode must be set and thatNode must NOT be set)
   @param insertionBranch -- the branch into which we insert   

 */
void insertNodeIntoBranch(state *chain, branch  toBeInserted, branch insertionBranch, double* blsNode1, double* blsNode2)
{  
  tree *tr = chain->tr; 
  int
    numBranches = getNumBranches( tr); 

  nodeptr
    pruned1 = findEmptyNodePtr(tr->nodep[toBeInserted.thisNode]); 

  assert(branchExists(tr, insertionBranch)); 
  
  nodeptr insert1 = findNodeFromBranch(tr, insertionBranch),
    insert2 = findNodeFromBranch(tr, invertBranch(insertionBranch)); 

  

  hookup(pruned1,insert1, blsNode1, numBranches );
  nodeptr pruned2 = findEmptyNodePtr( tr->nodep[toBeInserted.thisNode] ); 
  hookup(pruned2, insert2, blsNode2, numBranches);   


#ifdef DEBUG_SHOW_TOPO_CHANGES
  printInfo(chain, "inserted  %d into %d (bl=%f)\t%d (bl=%f)\n", pruned1->number, insert1->number,  insert2->number,pruned1->z[0], insert2->z[0]); 
#endif
}



/**
   @brief inverts the orientation of a branch 
 */ 
branch invertBranch(branch b)
{
  branch result ; 
  result.thisNode = b.thatNode; 
  result.thatNode = b.thisNode; 
  return result; 
}




/**
   @brief Gets a branch struct from two ids. 
 */ 
branch constructBranch(int thisNode, int thatNode)
{
  branch result; 
  result.thisNode = thisNode; 
  result.thatNode = thatNode; 

  return result; 
}




/**
   @brief Indicates  whether true nodes are hooked.  
 */
boolean branchExists(tree *tr, branch b)
{
  nodeptr a = tr->nodep[b.thisNode]; 

  return a->back->number == b.thatNode
    || a->next->back->number == b.thatNode
    || a->next->next->back->number == b.thatNode;   
}







/**
   @brief finds the correct nodeptr for a node that is currently
   hooked to another node.

   We should rather use node-id based stuff instead of pointers for
   hooking up things. This is more error-proof, since as soon as
   something is not hooked up as expected this function will trigger
   an assertion.
   
   @param targetNode   the node the nodeptr belongs to  
   @param neighborNode the id of the node, the target node should be hooked with
 */
static nodeptr findNodePtrByNeighbor(tree *tr, int targetNode,  int neighborNode)
{
  nodeptr 
    result = tr->nodep[targetNode];
  nodeptr
    p = result; 

  assert(targetNode < 2 * tr->mxtips); 
  assert(neighborNode < 2* tr->mxtips); 

  do 
    {
      if(p->back->number == neighborNode)
	{
	  assert(p->number == targetNode); 
	  return p; 
	}
      p = p->next; 
    } while(p != result) ; 

  /* this of course should never happen */
  assert(0); 
  return NULL;   
}


/* TODO replaces the next method   */
nodeptr findNodeFromBranch(tree *tr, branch b )
{  
  return findNodePtrByNeighbor(tr, b.thisNode, b.thatNode); 
}




/**
   @brief finds the root in the associated tree  
 */ 
branch findRoot(state *chain)
{
  tree *tr = chain->tr; 
  branch root = {0,0}; 
  for(int i = tr->mxtips +1 ; i < 2* tr->mxtips-1 ; ++i)
    {
      nodeptr
	p = tr->nodep[i],
	q = p;       
      do 
	{
	  if(q->x && q->back->x)
	    {
	      root.thisNode = q->number; 
	      root.thatNode = q->back->number; 	      
	    }
	  q = q->next; 
	} while(p != q); 
    }


  for(int i = 1; i < tr->mxtips+1; ++i)
    {
      nodeptr
	p = tr->nodep[i]; 
      if(p->back->x)
	{
	  root.thisNode = p->number; 
	  root.thatNode = p->back->number; 
	}
    }

  return root; 
}



/**
   @brief prunes a branch from the tree 
   @param  b -- the branch to be pruned (only thisNode is the pruning point, we only know thatNode for orientation)
   @param z -- the branch lengths at the pruning point after pruning 
*/
void pruneBranch(state *chain, branch b, double *z)
{  
  tree *tr = chain->tr ; 
  int numBl = getNumBranches( chain->tr); 

  nodeptr toPrune  = findNodeFromBranch(tr, b); 
  nodeptr pruned1 = toPrune->next->back,
    pruned2 = toPrune->next->next->back; 

  toPrune->next->back = toPrune->next->next->back = (nodeptr)NULL; 

  hookup(pruned1 , pruned2, z, numBl);

#ifdef DEBUG_SHOW_TOPO_CHANGES
  printInfo(chain, "pruning %d from %d,%d\tnew bl=%f\n", toPrune->number, pruned1->number, pruned2->number, z[0]);   
#endif
  
}
