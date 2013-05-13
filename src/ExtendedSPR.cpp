#include "ExtendedSPR.hpp"
#include "output.h"
#include "Path.hpp"
#include "TreeAln.hpp"
#include "eval.h"


// #define DEBUG_ESPR


ExtendedSPR::ExtendedSPR(double _relativeWeight, double _stopProb, double _multiplier)
  : stopProb(_stopProb), multiplier(_multiplier)    
{
  this->relativeProbability = _relativeWeight; 
  this->name = "eSPR"; 
  category = TOPOLOGY; 
  // ptype = E_SPR; 
}



ExtendedSPR::~ExtendedSPR()
{
}



// Meh 
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
void ExtendedSPR::destroyOrientationAlongPath( Path& path, tree *tr,  nodeptr p)
{  
  /* TODO efficiency =/  */

  if(NOT path.nodeIsOnPath(p->number) || isTip(p->number, tr->mxtips))
    return; 

  disorientHelper(tr,p);
  destroyOrientationAlongPath(path, tr, p->next->back); 
  destroyOrientationAlongPath(path, tr, p->next->next->back);
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

  assert(modifiedPath.size( ) == 0 ); 
  tree *tr  = traln.getTr();   

  branch start; 
  nodeptr p,q,r; 
  do 
    {
      start  = rand.drawSubtreeUniform(traln); 
  
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

  modifiedPath.append(start); 
  modifiedPath.append(constructBranch(q->number, r->number)); 

  nodeptr currentNode = rand.drawRandDouble01() ? q : r; 
  boolean accepted = FALSE;   
  while(NOT accepted)
    {
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      modifiedPath.pushToStackIfNovel(constructBranch(currentNode->number, n->number),tr->mxtips); 

      currentNode = n; 
      
      accepted = rand.drawRandDouble01() < stopProp && modifiedPath.size() > 2 ; 	
    }

  /* undo changes to the tree  */
  hookup(p->next,q, q->z, traln.getNumBranches()); 
  hookup(p->next->next,r,zqr, traln.getNumBranches());   

  /* now correct  */
  if(nodeIsInBranch(modifiedPath.at(1).thisNode, modifiedPath.at(2) ))    
    modifiedPath.at(1).thatNode = p->number; 
  else if(nodeIsInBranch(modifiedPath.at(1).thatNode, modifiedPath.at(2) ))    
    modifiedPath.at(1).thisNode = p->number; 
  else 
    assert(0); 

  /* correct the incorrectly set first branch in the path */
  modifiedPath.at(0) = getThirdBranch(tr, modifiedPath.at(0), modifiedPath.at(1)); 
  
#ifdef DEBUG_ESPR
  cout << *modifiedPath << endl; 
#endif
  
  modifiedPath.debug_assertPathExists(traln); 
}


/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
void ExtendedSPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  debug_printTree(traln);

  modifiedPath.clear(); 
  assert(modifiedPath.size() == 0); 
  drawPathForESPR( traln,rand ,stopProb); 

  modifiedPath.saveBranchLengthsPath(traln); 

  applyPathAsESPR(traln);

#ifdef ESPR_MULTIPLY_BL
  if(not prior.believingInFixedBranchLengths())
    {
      multiplyAlongBranchESPR(traln, rand, hastings, prior);
    }
#ifdef DEBUG_SHOW_TOPO_CHANGES
  cout << "after multiply: " << traln  << endl; 
#endif
#endif

  debug_checkTreeConsistency(traln.getTr()); 
}



void ExtendedSPR::resetState(TreeAln &traln, PriorBelief &prior )
{
  resetAlongPathForESPR (traln, prior);   

  /* TODO resetAlong... should be able to correctly restore branch lengths    */
  modifiedPath.restoreBranchLengthsPath(traln, prior ); 

  debug_checkTreeConsistency(traln.getTr());

  debug_printTree(traln); 
}


void ExtendedSPR::evaluateProposal(TreeAln &traln, PriorBelief &prior)
{  
  tree *tr = traln.getTr(); 

  branch futureRoot = getThirdBranch(tr, modifiedPath.at(0), modifiedPath.at(1)); 
  
  /* evaluate at root of inserted subtree */
  nodeptr toEval = findNodeFromBranch(tr, futureRoot); /* dangerous */

  destroyOrientationAlongPath(modifiedPath, tr, toEval); 
  destroyOrientationAlongPath(modifiedPath, tr, toEval->back);

  evaluateGenericWrapper(traln, toEval, FALSE);
}


/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void ExtendedSPR::multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior )
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
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void ExtendedSPR::applyPathAsESPR(TreeAln &traln )
{
  int numBranches = traln.getNumBranches(); 
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
  hookup(prPtr, prNPtr, prNPtr->z,numBranches); 
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL; 

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  nodeptr iPtr = findNodeFromBranch(tr, constructBranch(s2LastNode, lastNode)),
    iNPtr = findNodeFromBranch(tr, constructBranch(lastNode, s2LastNode)); 

  hookup(toBeInsertedPath, iPtr, iPtr->z, numBranches); 
  hookup(toBeInserted, iNPtr, toBeInserted->z,numBranches); 
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



/**
   @brief undoes topological changes, if an eSPR move was done
   according to rPath.

   Only resets the spr move, not any branch length changes due to BL
   multiplying.
 */ 
void ExtendedSPR::resetAlongPathForESPR(TreeAln &traln, PriorBelief &prior)
{
  // cout << "before reset: " <<  traln << endl; 
  // modifiedPath.printWithBLs(traln);
  // cout << endl; 

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

  hookup(lNPtr, oNPtr, oNPtr->z, numBranches);   
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

  hookup(sTPtr->next,        fPtr, fPtr->z, numBranches ); 
  hookup(sTPtr->next->next,  sPtr, &ztmp, numBranches ); 

#ifdef DEBUG_SHOW_TOPO_CHANGES
  printf("RESET: inserting %d into %d (%.2f),%d (%.2f)\n", sTPtr->number, fPtr->number, traln.getBranchLength( fPtr,0), sPtr->number, sPtr->z[0]); 
#endif
}


AbstractProposal* ExtendedSPR::clone() const
{
  // tout << "cloning "  << name << endl;
  return new ExtendedSPR(relativeProbability, stopProb, multiplier);
}
