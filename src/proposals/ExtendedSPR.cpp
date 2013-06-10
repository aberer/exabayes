#include "ExtendedSPR.hpp"
#include "output.h"
#include "Path.hpp"
#include "TreeAln.hpp"
#include "eval.h"

double ExtendedSPR::relativeWeight = 5.;

// #define DEBUG_ESPR

ExtendedSPR::ExtendedSPR( double _stopProb, double _multiplier)
  : stopProb(_stopProb), multiplier(_multiplier)    
{
  this->name = "eSPR"; 
  category = TOPOLOGY; 
}


ExtendedSPR::~ExtendedSPR()
{
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

  // cout << "subtree " << start  << endl; 

  /* save branches and prune */
  double zqr[NUM_BRANCHES]; 
  for(int i = 0; i < traln.getNumBranches(); ++i)
    zqr[i] = traln.getBranchLength( r,0);   
  traln.clipNode(q,r, q->z[0]);
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
  traln.clipNode(p->next,q, q->z[0]); 
  traln.clipNode(p->next->next,r,zqr[0]);   

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

  move.applyPathAsESPR(traln, modifiedPath);

  // cout << modifiedPath << endl; 

  assert(traln.getNumBranches() == 1 ); 

  bool modifiesBl = false; 
  for(auto v : secVar)
    modifiesBl |= v.getCategory() == BRANCH_LENGTHS; 
  if(modifiesBl)
    {
      auto brPr = secVar.at(0).getPrior();
      move.multiplyAlongBranchESPR(traln, rand, hastings, prior, modifiedPath, multiplier, brPr);
    }

  debug_checkTreeConsistency(traln.getTr()); 
}



void ExtendedSPR::resetState(TreeAln &traln, PriorBelief &prior )
{
  move.resetAlongPathForESPR (traln, prior, modifiedPath);   
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

  move.destroyOrientationAlongPath(modifiedPath, tr, toEval); 
  move.destroyOrientationAlongPath(modifiedPath, tr, toEval->back);

  evaluateGenericWrapper(traln, toEval, FALSE);
}


AbstractProposal* ExtendedSPR::clone() const
{
  return new ExtendedSPR( *this); 
}
