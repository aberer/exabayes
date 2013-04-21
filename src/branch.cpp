#include "axml.h"
#include "bayes.h"
#include "branch.h" 
#include "adapters.h"
#include "output.h"  
#include "TreeAln.hpp" 


/* #define GUIDE_SPR_BRANCH */



/**
   @brief divide a branch employing a ratio   
   
   notice: resultA is directed to thisNode in the branch 
 */ 
void divideBranchLengthsWithRatio(tree *tr, double orig,  double ratio, double *resultA, double *resultB)
{
  *resultA = pow(orig, (double)STRETCH_FACTOR * ratio) ; 
  *resultB = pow(orig, (double)STRETCH_FACTOR * (1-ratio)); 

#ifdef GUIDE_SPR_BRANCH
  printf("divided %g by %g into\t%g,%g\n",  orig, ratio, branchLengthToReal(tr, *resultA), branchLengthToReal(tr, *resultB));
#endif
}



/**
   @brief combine to branch lengths
 */ 
double combineBranchLengths(tree *tr, double origA, double origB)
{  
  double result = pow(origA * origB, 1. / (double) STRETCH_FACTOR );

#ifdef GUIDE_SPR_BRANCH
  printf("combined %g,%g into %g\n", branchLengthToReal(tr, origA), branchLengthToReal(tr, origB), branchLengthToReal(tr, result)); 
#endif
  return result; 
}


/**
    @brief gets the ratio of a to the sum of a and b 

    @notice direction is important!
 */ 
double getRatio(tree *tr, double a, double b )
{
  double ratio = log(a) / (log(a) + log(b)) ; 

  if(NOT (0 < ratio && ratio < 1.))
    {
      printf("\n\nWARNING: %g,%g (%g,%g) in ratio %g\n", a,b,branchLengthToReal(tr, a),branchLengthToReal(tr,b), ratio);
    }
    
  return ratio ; 

}





int getOtherNode(int node, branch b)
{
  if(b.thisNode == node)
    return b.thatNode; 
  else 
    {
      assert(node == b.thatNode); 
      return b.thisNode;
    }    
}



/**
   @brief converts the raxml bl into the real value for printing
 */ 
double branchLengthToReal(tree *tr, double internalBL)
{
  // TODO also own class for that 
  // assert(getNumBranches(tr) == 1 ); 
  return -log(internalBL) * tr->fracchange; 
}


/**
   @brief converts the real branch length to the raxml internal value  
 */ 
double branchLengthToInternal(tree *tr, double realBL)
{
  // TODO also own class for that 
  // assert(getNumBranches(tr) == 1 ); 
  return exp(-(realBL / tr->fracchange)); 
}



/**
   @brief gets equality of branches irregardless of their orientation
 */
boolean branchEqualUndirected(branch b1, branch b2)
{
  return (b1.thisNode == b2.thisNode && b1.thatNode == b2.thatNode) 
    || (b1.thisNode == b2.thatNode && b1.thatNode == b2.thisNode); 
}

boolean isTipBranch(branch b, int numTip)
{
  return b.thisNode <= numTip 
    || b.thatNode <= numTip; 
}

boolean branchesAreConnected(branch b1, branch b2)
{
  /* TODO  */
  assert(0); 
}


/**
   @brief Gets the third branch that is connected to a node.   
 */ 
branch getThirdBranch(tree *tr, branch b1, branch b2)
{  
  int node = getIntersectingNode(b1,b2); 
  
  assert(NOT isTip(node, tr->mxtips)); 

  nodeptr
    p = tr->nodep[node] ,
    q = p; 
  do 
    {
      if(NOT nodeIsInBranch(q->back->number,b1)
	 && NOT nodeIsInBranch(q->back->number,b2))
	return constructBranch(q->number, q->back->number); 
      q = q->next; 
    }while(p != q); 
  
  assert(0); 
  return constructBranch(0,0); 
}




int getIntersectingNode(branch b1, branch b2)
{
  if(b1.thisNode == b2.thisNode 
     || b1.thisNode == b2.thatNode )
    return b1.thisNode; 
  else if(b1.thatNode == b2.thisNode 
	  || b1.thatNode == b2.thatNode)
    return b1.thatNode; 
  else 
    {
      assert(0); 
      return 0 ; 
    }    
}


boolean nodeIsInBranch(int number, branch b )
{
  return number == b.thisNode || number == b.thatNode; 
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

   @param insertBranch -- the branch we insert (thisNode must be set and thatNode must NOT be set)
   @param insertionBranch -- the branch into which we insert   

 */
void insertNodeIntoBranch(state *chain, branch toBeInserted, branch insertionBranch, double* blsNode1, double* blsNode2)
{  
  tree *tr = chain->traln->getTr(); 
  int
    numBranches = chain->traln->getNumBranches(); 

  nodeptr
    pruned1 = findEmptyNodePtr(tr->nodep[toBeInserted.thisNode]); 

  assert(branchExists(tr, insertionBranch)); 
  
  nodeptr insert1 = findNodeFromBranch(tr, insertionBranch),
    insert2 = findNodeFromBranch(tr, invertBranch(insertionBranch)); 

  

  hookup(pruned1,insert1, blsNode1, numBranches );
  nodeptr pruned2 = findEmptyNodePtr( tr->nodep[toBeInserted.thisNode] ); 
  hookup(pruned2, insert2, blsNode2, numBranches);   


#ifdef DEBUG_SHOW_TOPO_CHANGES
  printInfo(chain, "inserted  %d into %d (bl=%f)\t%d (bl=%f)\n", pruned1->number, insert1->number,  insert2->number, traln->getBranchLength( pruned1->number,0), traln->getBranchLength( insert2->number,0)); 
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

  return (a->back != NULL &&  a->back->number == b.thatNode) 
    || (a->next->back != NULL && a->next->back->number == b.thatNode ) 
    || ( a->next->next->back != NULL && a->next->next->back->number == b.thatNode ) ;  
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
      if(p->back && p->back->number == neighborNode)
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
branch findRoot(tree *tr)
{
  /* tree *tr = chain->tr;  */
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
	  assert(root.thisNode == 0); 
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
  tree *tr = chain->traln->getTr() ; 
  int numBl = chain->traln->getNumBranches(); 

  nodeptr toPrune  = findNodeFromBranch(tr, b); 
  nodeptr pruned1 = toPrune->next->back,
    pruned2 = toPrune->next->next->back; 

  toPrune->next->back = toPrune->next->next->back = (nodeptr)NULL; 

  hookup(pruned1 , pruned2, z, numBl);

#ifdef DEBUG_SHOW_TOPO_CHANGES
  printInfo(chain, "pruning %d from %d,%d\tnew bl=%f\n", toPrune->number, pruned1->number, pruned2->number, z[0]);   
#endif
  
}




ostream& operator<<(ostream& rhs, const branch &b )
{
  rhs << b.thisNode << "/" << b.thatNode << "(" << b.length[0] << ")" ; 
  return rhs; 
}