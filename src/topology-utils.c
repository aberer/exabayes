/**
   @file  topology-utils.c
   
   @brief Various functions for manipulation of topologies and branch lengths. 
*/ 
#include "axml.h"
#include "proposalStructs.h"
#include "globals.h"
#include "randomness.h"
#include "main-common.h"
#include "adapters.h"
#include "output.h"


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
nodeptr findNodePtrByNeighbor(tree *tr, int targetNode,  int neighborNode)
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

   @param insertNode -- id of the node to be inserted. 
   
 */
void insertNodeIntoBranch(state *chain, int insertNode, int branchNode1, int branchNode2, double* blsNode1, double* blsNode2)
{  
  tree *tr = chain->tr; 
  int numBranches = getNumBranches( tr); 

  nodeptr
    pruned1 = findEmptyNodePtr(tr->nodep[insertNode]); 
  
  nodeptr
    insert1 = findNodePtrByNeighbor(tr, branchNode1, branchNode2),
    insert2 = findNodePtrByNeighbor(tr, branchNode2, branchNode1); 

  hookup(pruned1,insert1, blsNode1, numBranches   );

  nodeptr pruned2 = findEmptyNodePtr( tr->nodep[insertNode] ); 
  hookup(pruned2, insert2, blsNode2, numBranches);   


#ifdef DEBUG_SHOW_TOPO_CHANGES
  printInfo(chain, "inserted  %d into %d (bl=%f)\t%d (bl=%f)\n", insertNode, branchNode1, branchNode2,insert1->z[0], insert2->z[0]); 
#endif
}



/**
   @brief prunes a node from two nodes it is hooked up to. 
   
   @param prunedNode -- the id of the node that we prune 
   @param neighbor -- a neighbor of this the prunedNode 
   @param nextNeighbor -- another neighbor of the prunedNode
   @param bls -- branch lengths for the branch between nextNeighbor and neighbor 
   
 */
void pruneNodeFromNodes(state *chain , int prunedNode,  int neighbor, int nextNeighbor, double *bls)
{
  int numBl = getNumBranches( chain->tr); 

  nodeptr tPtr = findNodePtrByNeighbor(chain->tr, prunedNode, neighbor); 
  nodeptr nb = tPtr->back; 
  nodeptr tPtr2 = findNodePtrByNeighbor(chain->tr, prunedNode, nextNeighbor); 
  nodeptr nnb = tPtr2->back; 
  
  assert(nb->number == neighbor); 
  assert(nnb->number == nextNeighbor); 

  assert(tPtr2->number == prunedNode && tPtr->number == prunedNode); 
  
  assert(tPtr->next == tPtr2 || tPtr->next->next == tPtr2); 

  tPtr2->back = tPtr->back = (nodeptr)NULL; 

  hookup(nb , nnb, bls, numBl);

#ifdef DEBUG_SHOW_TOPO_CHANGES
  printInfo(chain, "pruning %d from %d,%d\tnew bl=%f\n", prunedNode, neighbor, nextNeighbor, bls[0]);   
#endif
}




/**
   @brief A check to ensure that the tree is still consistent. 
 */ 
void traverseAndCount(nodeptr p, int *count, tree *tr )
{
  nodeptr q;  

  *count += 1;

  if (! isTip(p->number,tr->mxtips)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  traverseAndCount(q->back, count, tr);
	  q = q->next;
	} 
    }
}




void insertWithGenericBL (nodeptr insertNode, nodeptr branchNode, double *insertZ, double *branchNodeZ, double *neighbourZ ,  int numBranches)
{
  nodeptr thirdNode = branchNode->back;   

  /* if(processID == 0) */
  /*   { */
  /*     printf("insertNode=%d, branchNode=%d, thirdNode=%d\n", insertNode->number, branchNode->number, thirdNode->number);  */
  /*   } */

  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
  hookup(insertNode, insertNode->back, insertZ, numBranches);
}

void insertWithUnifBL (state* chain, nodeptr insertNode, nodeptr branchNode, int numBranches)
{
  /* BUG: drawRandDouble does not take an argument  */
  assert(0); 

  double r;
  nodeptr thirdNode = branchNode->back;
  double branchNodeZ[numBranches];
  double neighbourZ[numBranches];
  
  for(int i=0; i<numBranches;i++)
  {
   /* r=drawRandDouble(branchNode->z[i]); */
   branchNodeZ[i]=r;
   neighbourZ[i]=1-r; 
  }
  
  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
}



void insertWithUnifBLScaled(nodeptr insertNode, nodeptr branchNode, double scale, int numBranches)
{
  /* BUG: drawRandDouble does not take an argument  */
  assert(0); 


  double r;
  nodeptr thirdNode = branchNode->back;
  double branchNodeZ[numBranches];
  double neighbourZ[numBranches];
  
  for(int i=0; i<numBranches;i++)
  {
   /* r=drawRandDouble(scale*branchNode->z[i]); */
   branchNodeZ[i]=r;
   neighbourZ[i]=1-r; 
  }
  
  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
}
