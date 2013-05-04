
#include "axml.h"
#include "bayes.h"		// 
#include "globals.h"
// #include "randomness.h"
// #include "main-common.h"
#include "adapters.h"
#include "output.h"
#include "path.h"


#include "TreeAln.hpp"

/**
   @brief gets the tree length. 

   p should be a tip. 
 */ 
double getTreeLength(TreeAln *traln, nodeptr p)
{
  tree *tr  = traln->getTr();
  if(isTip(p->number, tr->mxtips))
    return traln->getBranchLength( p,0); 
  else 
    {
      return traln->getBranchLength( p,0) * getTreeLength(traln, p->next->back)
	* getTreeLength(traln, p->next->next->back);       
    }
  
  assert(0);
  return 0; 
}

void restoreBranchLengthsPath(TreeAln *traln, path *s)
{
  int numBranches = traln->getNumBranches();
  for(int i = 0; i < s->index; ++i)
    {
      branch b = s->content[i]; 
      nodeptr p = findNodeFromBranch(traln->getTr(), b); 
      hookup(p, p->back, b.length, numBranches); 
    }

}




static void disorientHelper(tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips))
    {

      /* printf("not disorienting tip %d\n", p->number);  */
    }
  else if(p->x)
    {
      p->x = 0; 
      p->next->x = 1; 
      /* printf("disorienting %d (CORRECT  before) -> oriented to %d NOW \n", p->number, p->next->back->number); */
    }
  else 
    {
      /* printf("disorienting %d  (was incorrect before)\n", p->number); */
    }
}



/**
   @brief dis-orients the path, s.t. the lnl can be recomputed
   correctly.
   
   @notice assumes that the lnl will be evaluated at the end of the
   path; also notice that an SPR move has already been applied to the
   tree
 */ 
void destroyOrientationAlongPath(tree *tr, path *rPath, nodeptr p)
{  
  /* TODO efficiency =/  */


  if(NOT nodeIsOnPath(p->number, rPath) || isTip(p->number, tr->mxtips))
    return; 

  disorientHelper(tr,p);
  destroyOrientationAlongPath(tr, rPath, p->next->back); 
  destroyOrientationAlongPath(tr, rPath, p->next->next->back);
}


/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void applyPathAsESPR(TreeAln *traln, path *rPath )
{
  int numBranches = traln->getNumBranches(); 
  tree *tr = traln->getTr();
  
  assert(stackLength(rPath) > 2 ); 

  /* get the subtree ptr */
  /* nodeptr sTPtr = findNodeFromBranch(tr, rPath->content[0]);  */
  nodeptr sTPtr = findNodeFromBranch(tr,getThirdBranch(tr, rPath->content[0], rPath->content[1])); 
  /* printf("sTPtr = %d,%d\n", sTPtr->number, sTPtr->back->number );  */
  assert(nodeIsInBranch(sTPtr->number, rPath->content[1])); 

  /* find the two pointers to hook up  */
  nodeptr prPtr = sTPtr->next->back,
    prNPtr = sTPtr->next->next->back; 
    

  /* swap, s.t. prPtr is the one that is oriented towards the path */
  if(nodeIsInBranch(prNPtr->number, rPath->content[1]))
    {
      nodeptr tmp = prPtr; 
      prPtr = prNPtr; 
      prNPtr = tmp; 
    }
  nodeptr toBeInserted = prPtr->back;   
  nodeptr toBeInsertedPath = prNPtr->back;  
  assert(toBeInserted->number == sTPtr->number && sTPtr->number == toBeInsertedPath->number); 

  
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("APPLY: pruning %d from %d,%d\n", toBeInserted->number, prPtr->number, prNPtr->number); 
#endif
  /* prune sTPtr */
  hookup(prPtr, prNPtr, prNPtr->z,numBranches); 
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL; 

  /* find insertion nodes */
  nodeptr iPtr = findNodeFromBranch(tr, rPath->content[rPath->index-1]),
    iNPtr = findNodeFromBranch(tr, invertBranch(rPath->content[rPath->index-1]));

  /* swap, if iPtr is not directed towards the path */
  if(nodeIsInBranch(iNPtr->number, rPath->content[rPath->index-1]))
    {
      nodeptr tmp = iNPtr; 
      iNPtr = iPtr; 
      iPtr = tmp; 
    }

  hookup(toBeInsertedPath, iPtr, iPtr->z, numBranches); 
  hookup(toBeInserted, iNPtr, toBeInserted->z,numBranches); 
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("APPLY: inserting %d to %d,%d\n", toBeInserted->number, iPtr->number, iNPtr->number); 
#endif
  
  assert(branchExists(tr, constructBranch(sTPtr->number, iNPtr->number))); 
  assert(branchExists(tr, constructBranch(sTPtr->number, iPtr->number))); 
}


/**
   @brief undoes topological changes, if an eSPR move was done
   according to rPath.

   Only resets the spr move, not any branch length changes due to BL
   multiplying.
 */ 
void resetAlongPathForESPR(TreeAln *traln, path *rPath)
{
  tree *tr = traln->getTr();

  /* BUG branch lengths are not restored correctly, currently relying on restoreBranchLengthsPath function */

  assert(rPath->index > 2); 
  int numBranches = traln->getNumBranches(); 
  int lastNode = 0, otherNode = 0; 
  if(nodeIsInBranch(rPath->content[rPath->index-1].thisNode,  rPath->content[rPath->index-2]))
    {
      lastNode = rPath->content[rPath->index-1].thatNode; 
      otherNode = rPath->content[rPath->index-1].thisNode; 
    }
  else 
    {
      lastNode = rPath->content[rPath->index-1].thisNode;  
      otherNode = rPath->content[rPath->index-1].thatNode; 
    }

  /* prune the subtree  */
  int sTNode = getIntersectingNode(rPath->content[0], rPath->content[1]); 
  nodeptr lNPtr = findNodeFromBranch(tr, constructBranch(lastNode, sTNode)),
    oNPtr = findNodeFromBranch(tr, constructBranch(otherNode, sTNode)); 
  
  double ztmp[NUM_BRANCHES]; 
  for(int i = 0; i < numBranches; ++i)
    ztmp[i] = traln->getBranchLength( lNPtr,0); 

  hookup(lNPtr, oNPtr, oNPtr->z, numBranches);   
  branch b = getThirdBranch(tr, constructBranch(sTNode, lastNode), constructBranch(sTNode, otherNode)); 
  assert(b.thisNode == sTNode); 

  nodeptr sTPtr = findNodeFromBranch(tr, b);
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL;
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("RESET: pruning %d from %d (%.2f),%d (%.2f)\n", sTPtr->number, oNPtr->number, traln->getBranchLength( oNPtr->number,0) , lNPtr->number, traln->getBranchLength( lNPtr->number,0)); 
#endif

  int firstNode = rPath->content[0].thisNode == sTNode ? rPath->content[0].thatNode : rPath->content[0].thisNode; 
  int secondNode = rPath->content[1].thisNode == sTNode ? rPath->content[1].thatNode : rPath->content[1].thisNode; 

  nodeptr fPtr = findNodeFromBranch(tr, constructBranch(firstNode,secondNode )),
    sPtr = findNodeFromBranch(tr, constructBranch(secondNode,firstNode )); 

  
  
  hookup(sTPtr->next,        fPtr, fPtr->z, numBranches ); 
  hookup(sTPtr->next->next,  sPtr, ztmp, numBranches ); 
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("RESET: inserting %d into %d (%.2f),%d (%.2f)\n", sTPtr->number, fPtr->number, traln->getBranchLength( fPtr->number,0), sPtr->number, sPtr->z[0]); 
#endif
  
}





/**
   @brief A check to ensure that the tree is still consistent. 
   
   @param  p -- pointer to an inner (!) node
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

  hookup(insertNode->next , branchNode, branchNodeZ, numBranches); 
  hookup(insertNode->next->next, thirdNode, neighbourZ, numBranches); 
  hookup(insertNode, insertNode->back, insertZ, numBranches);
}

void insertWithUnifBL (Chain* chain, nodeptr insertNode, nodeptr branchNode, int numBranches)
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

