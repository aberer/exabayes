#include "SprMove.hpp"

// Meh 
static void disorientHelper(tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips))
    {
      // printf("not disorienting tip %d\n", p->number);
    }
  else if(p->x)
    {
      p->x = 0; 
      p->next->x = 1; 
      // printf("disorienting %d (CORRECT  before) -> oriented to %d NOW \n", p->number, p->next->back->number);
    }
  else 
    {
      // printf("disorienting %d  (was incorrect before)\n", p->number);
    }
}



/**
   @brief dis-orients the path, s.t. the lnl can be recomputed
   correctly.
   
   @notice assumes that the lnl will be evaluated at the end of the
   path; also notice that an SPR move has already been applied to the
   tree
 */ 
void SprMove::destroyOrientationAlongPath( Path& path, tree *tr,  nodeptr p)
{  
  /* TODO efficiency =/  */

  if(NOT path.nodeIsOnPath(p->number) || isTip(p->number, tr->mxtips))
    return; 

  disorientHelper(tr,p);
  destroyOrientationAlongPath(path, tr, p->next->back); 
  destroyOrientationAlongPath(path, tr, p->next->next->back);
}



/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void SprMove::multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior,  Path &modifiedPath, double multiplier )
{
  assert(modifiedPath.size() >= 2); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
 
  int sTNode = modifiedPath.getNthNodeInPath(1); 
  branch firstBranch = constructBranch( modifiedPath.getNthNodeInPath(0), modifiedPath.getNthNodeInPath(2)); 

  modifiedPath.multiplyBranch(traln, rand, firstBranch, multiplier, hastings, prior); 
  
  /* treat all branches except the first 2 and the last one */
  int ctr = 0; 
  for(int i = 0; i < modifiedPath.size() ; ++i)
    {
      if(ctr < 2 )
	continue; 
      
      branch &b = modifiedPath.at(i); 

      modifiedPath.multiplyBranch(traln, rand, b, multiplier, hastings, prior); 
    }

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  modifiedPath.multiplyBranch(traln, rand, constructBranch(sTNode, lastNode), multiplier, hastings, prior); 
  modifiedPath.multiplyBranch(traln, rand, constructBranch(sTNode, s2LastNode), multiplier, hastings, prior); 
}



/**
   @brief undoes topological changes, if an eSPR move was done
   according to rPath.

   Only resets the spr move, not any branch length changes due to BL
   multiplying.
 */ 
void SprMove::resetAlongPathForESPR(TreeAln &traln, PriorBelief &prior, Path& modifiedPath)
{
  tree *tr = traln.getTr();

  /* BUG branch lengths are not restored correctly, currently relying on restoreBranchLengthsPath function */

  assert(modifiedPath.size() > 2); 
  int numBranches = traln.getNumBranches(); 

  int numNodes = modifiedPath.getNumberOfNodes(); 
  int lastNode = modifiedPath.getNthNodeInPath(numNodes-1),
    otherNode = modifiedPath.getNthNodeInPath(numNodes-2); 

  /* prune the subtree  */  
  int sTNode =  modifiedPath.getNthNodeInPath(1); 
  nodeptr lNPtr = findNodeFromBranch(tr, constructBranch(lastNode, sTNode)),
    oNPtr = findNodeFromBranch(tr, constructBranch(otherNode, sTNode)); 
  
  // double ztmp[NUM_BRANCHES]; 
  double ztmp ; 
  for(int i = 0; i < numBranches; ++i)
    ztmp = traln.getBranchLength( lNPtr,0); 

  traln.clipNode(lNPtr, oNPtr, oNPtr->z[0]);   
  branch b = getThirdBranch(tr, constructBranch(sTNode, lastNode), constructBranch(sTNode, otherNode)); 
  assert(b.thisNode == sTNode); 

  nodeptr sTPtr = findNodeFromBranch(tr, b);
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL;
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("RESET: pruning %d from %d (%.2f),%d (%.2f)\n", sTPtr->number, oNPtr->number, traln.getBranchLength( oNPtr,0) , lNPtr->number, traln.getBranchLength( lNPtr,0)); 
#endif

  int firstNode = modifiedPath.getNthNodeInPath(0),
    secondNode = modifiedPath.getNthNodeInPath(2); 

  nodeptr fPtr = findNodeFromBranch(tr, constructBranch(firstNode,secondNode )),
    sPtr = findNodeFromBranch(tr, constructBranch(secondNode,firstNode )); 

  traln.clipNode(sTPtr->next,        fPtr, fPtr->z[0] ); 
  traln.clipNode(sTPtr->next->next,  sPtr, ztmp ); // 

#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("RESET: inserting %d into %d (%.2f),%d (%.2f)\n", sTPtr->number, fPtr->number, traln.getBranchLength( fPtr,0), sPtr->number, sPtr->z[0]); 
#endif
}



/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void SprMove::applyPathAsESPR(TreeAln &traln, Path &modifiedPath )
{
  tree *tr = traln.getTr();

#ifdef CONTROL_ESPR
double treeLengthBefore = traln.getTreeLength(); 
#endif
  
  assert(modifiedPath.size() > 2 ); 

  /* get the subtree ptr */
  nodeptr sTPtr = findNodeFromBranch(tr,getThirdBranch(tr, modifiedPath.at(0), modifiedPath.at(1))); 
  assert(nodeIsInBranch(sTPtr->number, modifiedPath.at(1))); 


  // finds the two nodeptrs adjacent to the subtree  
  nodeptr prNPtr = findNodeFromBranch(tr, constructBranch(modifiedPath.getNthNodeInPath(0) ,modifiedPath.getNthNodeInPath(1))),
    prPtr = findNodeFromBranch(tr, constructBranch(modifiedPath.getNthNodeInPath(2), modifiedPath.getNthNodeInPath(1))); 

  nodeptr toBeInserted = prPtr->back;   
  nodeptr toBeInsertedPath = prNPtr->back;  
  assert(toBeInserted->number == sTPtr->number && sTPtr->number == toBeInsertedPath->number); 

  
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("APPLY: pruning %d from %d,%d\n", toBeInserted->number, prPtr->number, prNPtr->number); 
#endif
  /* prune sTPtr */
  traln.clipNode(prPtr, prNPtr, prNPtr->z[0]); 
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL; 

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  nodeptr iPtr = findNodeFromBranch(tr, constructBranch(s2LastNode, lastNode)),
    iNPtr = findNodeFromBranch(tr, constructBranch(lastNode, s2LastNode)); 

  traln.clipNode(toBeInsertedPath, iPtr, iPtr->z[0]); 
  traln.clipNode(toBeInserted, iNPtr, toBeInserted->z[0]); 
#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("APPLY: inserting %d to %d,%d\n", toBeInserted->number, iPtr->number, iNPtr->number); 
#endif
  
  assert(branchExists(tr, constructBranch(sTPtr->number, iNPtr->number))); 
  assert(branchExists(tr, constructBranch(sTPtr->number, iPtr->number))); 

#ifdef DEBUG_SHOW_TOPO_CHANGES
  cout << "after inserting: " << *traln<< endl; 
#endif

#ifdef CONTROL_ESPR
  double treeLengthAfter =  traln.getTreeLength(); 
  if( fabs(treeLengthAfter  -  treeLengthBefore) > 1e-3  )
    {
      cout << setprecision(8)  << "TL before " << branchLengthToReal(traln.getTr(),treeLengthBefore) << "\tafter" <<  branchLengthToReal(traln.getTr(), treeLengthAfter) << endl; 
      assert(treeLengthAfter == treeLengthBefore); 
    }
#endif
}
