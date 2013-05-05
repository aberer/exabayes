#include "ExtendedSPR.hpp"
#include "output.h"
#include "Path.hpp"
#include "TreeAln.hpp"


ExtendedSPR::ExtendedSPR(Chain *_chain, double _relativeWeight, double _stopProb, double _multiplier)
  : chain(_chain), stopProb(_stopProb), multiplier(_multiplier)    
{
  this->relativeProbability = _relativeWeight; 
  this->name = "eSPR"; 
  category = TOPOLOGY; 
  ptype = E_SPR; 
  modifiedPath = new Path(); 
}



ExtendedSPR::~ExtendedSPR()
{
  delete modifiedPath; 
}


/**
   @brief draws a random set of branches in the tree that constitute a
   path
   
   This function employs the eSPR strategy 
   
   @param s -- the result : the first two branches in the stack define
   the root of the pruned subtree
 */
void ExtendedSPR::drawPathForESPR(TreeAln& traln, Randomness &rand, double stopProp )
{
  assert(0 < stopProp && stopProp < 1.0 ); 

  assert(modifiedPath->stackIsEmpty()); 
  tree *tr  = traln.getTr();   

  branch start; 
  nodeptr p,q,r; 
  do 
    {
      start  = rand.drawSubtreeUniform(&traln); 
  
      p = findNodeFromBranch(tr, start);
      q = p->next->back; 
      r = p->next->next->back;
    } while(isTip(q->number,tr->mxtips) && isTip(r->number, tr->mxtips)); 

  /* TODO assert that the pruned subtree is not too big (s.t. we cannot find any insertion points)  */


  /* save branches and prune */
  double zqr[NUM_BRANCHES]; 
  for(int i = 0; i < traln.getNumBranches(); ++i)
    zqr[i] = traln.getBranchLength( r,0); 
  hookup(q,r, q->z, traln.getNumBranches());
  p->next->back = p->next->next->back = (nodeptr)NULL; 

  modifiedPath->pushStack(start); 
  modifiedPath->pushStack(constructBranch(q->number, r->number)); 

  nodeptr currentNode = rand.drawRandDouble01() ? q : r; 
  boolean accepted = FALSE;   
  while(NOT accepted)
    {
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      modifiedPath->pushToStackIfNovel(constructBranch(currentNode->number, n->number),tr->mxtips); 

      currentNode = n; 
      
      accepted = rand.drawRandDouble01() < stopProp && modifiedPath->stackLength() > 2 ; 	
#ifdef DEBUG_ESPR
      cout << *modifiedPath << endl; 
#endif
    }

  /* undo changes to the tree  */
  hookup(p->next,q, q->z, traln.getNumBranches()); 
  hookup(p->next->next,r,zqr, traln.getNumBranches());   

  /* now correct  */
  if(nodeIsInBranch(modifiedPath->at(1).thisNode, modifiedPath->at(2) ))    
    modifiedPath->at(1).thatNode = p->number; 
  else if(nodeIsInBranch(modifiedPath->at(1).thatNode, modifiedPath->at(2) ))    
    modifiedPath->at(1).thisNode = p->number; 
  else 
    assert(0); 

  /* correct the incorrectly set first branch in the path */
  modifiedPath->at(0) = getThirdBranch(tr, modifiedPath->at(0), modifiedPath->at(1)); 
  
#ifdef DEBUG_ESPR
  cout << *modifiedPath << endl; 
#endif
  
  modifiedPath->debug_assertPathExists(traln); 
}


/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
void ExtendedSPR::applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand)
{
  debug_printTree(chain);

  modifiedPath->clearStack(); 
  assert(modifiedPath->stackLength() == 0); 
  drawPathForESPR( traln,rand ,stopProb); 

  modifiedPath->saveBranchLengthsPath(traln); 

  applyPathAsESPR(chain->traln);

#ifdef ESPR_MULTIPLY_BL
  multiplyAlongBranchESPR(traln, rand, multiplier, hastings);
#endif

  debug_checkTreeConsistency(chain->traln->getTr()); 
}



void ExtendedSPR::resetState(TreeAln &traln, PriorManager &prior )
{
  resetAlongPathForESPR (&traln);   

  /* TODO resetAlong... should be able to correctly restore branch lengths    */
  modifiedPath->restoreBranchLengthsPath(traln); 

  debug_checkTreeConsistency((&traln)->getTr());

  debug_printTree(chain); 
}




void ExtendedSPR::evaluateProposal(TreeAln &traln, PriorManager &prior)
{  
  tree *tr = traln.getTr(); 
  Path *rPath = modifiedPath; 

  branch futureRoot = getThirdBranch(tr, modifiedPath->at(0), modifiedPath->at(1)); 
  
  /* evaluate at root of inserted subtree */
  nodeptr toEval = findNodeFromBranch(tr, futureRoot); /* dangerous */

  rPath->destroyOrientationAlongPath(tr, toEval); 
  rPath->destroyOrientationAlongPath(tr, toEval->back);

  /* printf("evaluating at branch %d,%d\n", toEval->number, toEval->back->number);  */
  evaluateGenericWrapper(chain, toEval, FALSE);
}





/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void ExtendedSPR::multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double multi, double &hastings )
{
  tree *tr = traln.getTr();

  assert(modifiedPath->stackLength() >= 2); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 

  /* first two branches: consider that subtree has been pruned */
  branch b = getThirdBranch(tr, modifiedPath->at(0), modifiedPath->at(1) ); 
  int sTNode = b.thisNode; 

  branch firstBranch = constructBranch( getOtherNode(sTNode, modifiedPath->at(0)) ,
					getOtherNode(sTNode, modifiedPath->at(1)));   

  modifiedPath->multiplyBranch(traln, rand, firstBranch, multi, &hastings); 
  
  /* treat all branches except the first 2 and the last one */
  int ctr = 0; 
  for(auto b : modifiedPath->getStack()  )
    {
      if(ctr < 2 )
	continue; 

      modifiedPath->multiplyBranch(traln, rand, b, multi, &hastings); 
    }

  /* treat last two paths (only one in representation) */
  branch lastB = modifiedPath->at(modifiedPath->stackLength()-1); 
  b = constructBranch(sTNode, lastB.thisNode);	// s->content[s->index-1]
  modifiedPath->multiplyBranch(traln, rand, b,multi, &hastings); 

  b = constructBranch(sTNode, lastB.thatNode); 
  modifiedPath->multiplyBranch(traln, rand,b,multi, &hastings);   
}


/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void ExtendedSPR::applyPathAsESPR(TreeAln *traln )
{
  int numBranches = traln->getNumBranches(); 
  tree *tr = traln->getTr();
  
  assert(modifiedPath->stackLength() > 2 ); 

  /* get the subtree ptr */
  nodeptr sTPtr = findNodeFromBranch(tr,getThirdBranch(tr, modifiedPath->at(0), modifiedPath->at(1))); 
  assert(nodeIsInBranch(sTPtr->number, modifiedPath->at(1))); 

  /* find the two pointers to hook up  */
  nodeptr prPtr = sTPtr->next->back,
    prNPtr = sTPtr->next->next->back; 
    

  /* swap, s.t. prPtr is the one that is oriented towards the path */
  if(nodeIsInBranch(prNPtr->number, modifiedPath->at(1)))
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
  nodeptr iPtr = findNodeFromBranch(tr, modifiedPath->at(modifiedPath->stackLength()-1)),
    iNPtr = findNodeFromBranch(tr, invertBranch(modifiedPath->at(modifiedPath->stackLength()-1)));

  /* swap, if iPtr is not directed towards the path */
  if(nodeIsInBranch(iNPtr->number, modifiedPath->at(modifiedPath->stackLength()-1)))
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
void ExtendedSPR::resetAlongPathForESPR(TreeAln *traln)
{
  tree *tr = traln->getTr();

  /* BUG branch lengths are not restored correctly, currently relying on restoreBranchLengthsPath function */

  assert(modifiedPath->stackLength() > 2); 
  int numBranches = traln->getNumBranches(); 
  // int lastNode = 0, otherNode = 0; 
  // branch &lastBranch = modifiedPath->at(modifiedPath->stackLength()-1),
  //   &secondToLastBranch = modifiedPath->at(modifiedPath->stackLength()-2); 

  int numNodes = modifiedPath->getNumberOfNodes(); 
  int lastNode = modifiedPath->getNthNodeInPath(numNodes),
    otherNode = modifiedPath->getNthNodeInPath(numNodes-1); 
  
  
  // if(nodeIsInBranch(lastBranch.thisNode, secondToLastBranch))
  //   {
  //     lastNode = lastBranch.thatNode; 
  //     otherNode = lastBranch.thisNode; 
  //   }
  // else 
  //   {
  //     lastNode = lastBranch.thisNode;  
  //     otherNode =lastBranch.thatNode; 
  //   }

  /* prune the subtree  */
  int sTNode = getIntersectingNode(modifiedPath->at(0), modifiedPath->at(1)); 
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
  printf("RESET: pruning %d from %d (%.2f),%d (%.2f)\n", sTPtr->number, oNPtr->number, traln->getBranchLength( oNPtr,0) , lNPtr->number, traln->getBranchLength( lNPtr,0)); 
#endif

  int firstNode = modifiedPath->at(0).thisNode == sTNode ? modifiedPath->at(0).thatNode : modifiedPath->at(0).thisNode; 
  int secondNode = modifiedPath->at(1).thisNode == sTNode ? modifiedPath->at(1).thatNode : modifiedPath->at(1).thisNode; 

  nodeptr fPtr = findNodeFromBranch(tr, constructBranch(firstNode,secondNode )),
    sPtr = findNodeFromBranch(tr, constructBranch(secondNode,firstNode )); 

  hookup(sTPtr->next,        fPtr, fPtr->z, numBranches ); 
  hookup(sTPtr->next->next,  sPtr, ztmp, numBranches ); 

#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("RESET: inserting %d into %d (%.2f),%d (%.2f)\n", sTPtr->number, fPtr->number, traln->getBranchLength( fPtr,0), sPtr->number, sPtr->z[0]); 
#endif
}
