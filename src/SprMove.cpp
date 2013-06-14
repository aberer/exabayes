#include <set>

#include "SprMove.hpp"


// #define CONTROL_ESPR

// Meh 
static void disorientHelper(const TreeAln &traln, nodeptr p)
{
  if(traln.isTipNode(p))
    {
    }
  else if(p->x)
    {
      p->x = 0; 
      p->next->x = 1; 
    }
  else 
    {
    }
}



/**
   @brief dis-orients the path, s.t. the lnl can be recomputed
   correctly.
   
   @notice assumes that the lnl will be evaluated at the end of the
   path; also notice that an SPR move has already been applied to the
   tree
 */ 
void SprMove::destroyOrientationAlongPath( Path& path, const TreeAln &traln,  nodeptr p)
{  
  nat first = path.getNthNodeInPath(0),
    last = path.getNthNodeInPath(path.getNumberOfNodes()-1); 

  if(NOT path.nodeIsOnPath(p->number) || traln.isTipNode(p) || p->number == first || p->number == last )
    return; 

  disorientHelper(traln, p);
  destroyOrientationAlongPath(path, traln, p->next->back); 
  destroyOrientationAlongPath(path, traln, p->next->next->back);
}



/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void SprMove::multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior,  Path &modifiedPath, double multiplier, shared_ptr<AbstractPrior> brPr)
{
  assert(modifiedPath.size() >= 2); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
 
  int sTNode = modifiedPath.getNthNodeInPath(1); 
  branch firstBranch = constructBranch( modifiedPath.getNthNodeInPath(0), modifiedPath.getNthNodeInPath(2)); 

  modifiedPath.multiplyBranch(traln, rand, firstBranch, multiplier, hastings, prior, brPr); 
  
  /* treat all branches except the first 2 and the last one */
  int ctr = 0; 
  for(nat i = 0; i < modifiedPath.size() ; ++i)
    {
      if(ctr < 2 )
	continue; 
      
      branch &b = modifiedPath.at(i); 
      modifiedPath.multiplyBranch(traln, rand, b, multiplier, hastings, prior, brPr); 
    }

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  modifiedPath.multiplyBranch(traln, rand, constructBranch(sTNode, lastNode), multiplier, hastings, prior, brPr); 
  modifiedPath.multiplyBranch(traln, rand, constructBranch(sTNode, s2LastNode), multiplier, hastings, prior, brPr); 
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

  assert(modifiedPath.size() > 2); 
  int numBranches = traln.getNumBranches(); 

  int numNodes = modifiedPath.getNumberOfNodes(); 
  int lastNode = modifiedPath.getNthNodeInPath(numNodes-1),
    otherNode = modifiedPath.getNthNodeInPath(numNodes-2); 

  /* prune the subtree  */  
  int sTNode =  modifiedPath.getNthNodeInPath(1); 
  nodeptr lNPtr = findNodeFromBranch(tr, constructBranch(lastNode, sTNode)),
    oNPtr = findNodeFromBranch(tr, constructBranch(otherNode, sTNode)); 

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
  if( fabs(treeLengthAfter  -  treeLengthBefore) > 1e-6  )
    {
      cout << setprecision(8)  << "TL before " << branchLengthToReal(traln.getTr(),treeLengthBefore) << "\tafter" <<  branchLengthToReal(traln.getTr(), treeLengthAfter) << endl; 
      assert(treeLengthAfter == treeLengthBefore); 
    }
#endif
}



/** 
    @brief Gets the description path after the move has been executed. 

    Can be used for reversal => orientation is inverted 

    Does not contain branch lengths.   
 */ 
void SprMove::getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath)
{  
  int subTreeNode =  modifiedPath.getNthNodeInPath(1);   
  branch lastBranch = modifiedPath.at(modifiedPath.size()-1); 

  resultPath.append(constructBranch( lastBranch.thatNode , subTreeNode ) )  ; 
  resultPath.append(constructBranch( lastBranch.thisNode , subTreeNode ) )  ; 
  
  assert(modifiedPath.size() > 2); 
  
  // insert the non-descriptive branches in reversed order 
  for(int i = int(modifiedPath.size())-2 ; i > 1 ; --i)
    resultPath.append(modifiedPath.at(i)); 

  // fuse the first two branches 
  auto b1 = modifiedPath.at(0),
    b2  = modifiedPath.at(1); 

  set<int> ids; 
  ids.insert(b1.thisNode); 
  ids.insert(b2.thisNode); 
  ids.insert(b1.thatNode); 
  ids.insert(b2.thatNode); 
  
  assert(ids.find(subTreeNode) != ids.end());   
  ids.erase(subTreeNode); 
  
  assert(ids.size() == 2); 

  
  auto it = ids.begin(); 
  int a = *it; 
  ++it; 
  int b = *it;   
  
  // auto newBranch = constructBranch(a,b); 
  // if(branchesAreConnected(newBranch,  resultPath.at(resultPath.size()-1)

  resultPath.append(constructBranch(a,b)); 

  assert(modifiedPath.size() == resultPath.size()); 
}
