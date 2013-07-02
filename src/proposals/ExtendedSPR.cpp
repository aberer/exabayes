#include "ExtendedSPR.hpp"
#include "Path.hpp"
#include "TreeAln.hpp"


// #define DEBUG_ESPR

ExtendedSPR::ExtendedSPR( double _stopProb, double _multiplier)
  : stopProb(_stopProb), multiplier(_multiplier)    
{
  this->name = "eSPR"; 
  category = Category::TOPOLOGY; 
  relativeWeight = 5.;
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
  Path modifiedPath; 
  assert(0 < stopProp && stopProp < 1.0 ); 

  assert(modifiedPath.size( ) == 0 ); 
  tree *tr  = traln.getTr();   

  Branch start; 
  nodeptr p,q,r; 
  do 
    {
      start = traln.drawBranchWithInnerNode(rand); 

      p = start.findNodePtr(traln );
      q = p->next->back; 
      r = p->next->next->back;
    } while(traln.isTipNode(q) || traln.isTipNode(r) ); 

  /* save branches and prune */
  double zqr[NUM_BRANCHES]; 
  for(int i = 0; i < traln.getNumBranches(); ++i)
    zqr[i] = traln.getBranchLength( r,0);   
  traln.clipNode(q,r, q->z[0]);
  p->next->back = p->next->next->back = (nodeptr)NULL; 

  modifiedPath.append(start); 
  modifiedPath.append(Branch(q->number, r->number)); 

  nodeptr currentNode = rand.drawRandDouble01() ? q : r; 
  boolean accepted = FALSE;   
  while(NOT accepted)
    {      
      nodeptr n = 
	rand.drawRandDouble01()  < 0.5 
	? currentNode->next->back
	: currentNode->next->next->back; 

      modifiedPath.pushToStackIfNovel(Branch(currentNode->number, n->number),traln); 

      currentNode = n; 
      
      accepted = rand.drawRandDouble01() < stopProp && modifiedPath.size() > 2 ; 	
    }

  /* undo changes to the tree  */
  traln.clipNode(p->next,q, q->z[0]); 
  traln.clipNode(p->next->next,r,zqr[0]);   

  /* now correct  */
  if( modifiedPath.at(2).nodeIsInBranch(modifiedPath.at(1).getPrimNode()))    
    modifiedPath.at(1).setSecNode(p->number); 
  else if(modifiedPath.at(2).nodeIsInBranch(modifiedPath.at(1).getSecNode()  ))    
    modifiedPath.at(1).setPrimNode( p->number); 
  else 
    assert(0); 

  /* correct the incorrectly set first branch in the path */
  modifiedPath.at(0) = modifiedPath.at(0).getThirdBranch(traln, modifiedPath.at(1)); 
  
#ifdef DEBUG_ESPR
  cout << *modifiedPath << endl; 
#endif
  
  modifiedPath.debug_assertPathExists(traln); 

  // TODO in principle, we can throw away the path of this proposal and use the move proposal  
  // move.

  // BEGIN  TODO remove modifiedPath
  Branch bla = modifiedPath.at(modifiedPath.size()-1); 
  move.extractMoveInfo(traln, {Branch(p->number, p->back->number), Branch(bla.getPrimNode(), bla.getSecNode()) }); 
  modifiedPath.clear(); 
  // END
}



/**
   @brief applies an extended spr move (our version)
   
   the same function as below, but I cleaned the other flavour, since so much stuff has changed. 
 */ 
void ExtendedSPR::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  // debug_printTree(traln);
  drawPathForESPR( traln,rand ,stopProb); 

  move.applyToTree(traln); 

  assert(traln.getNumBranches() == 1 ); 

  bool modifiesBl = false; 
  for(auto v : secVar)
    modifiesBl |= v->getCategory() == Category::BRANCH_LENGTHS; 

#ifdef NO_SEC_BL_MULTI
  modifiesBl = false; 
#endif

  if(modifiesBl)
    {
      auto brPr = secVar.at(0)->getPrior();
      move.multiplyBranches(traln, rand, hastings, prior, multiplier, {brPr} ); 
    }

  // debug_checkTreeConsistency(traln); 
}


void ExtendedSPR::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior)
{  
  Branch toEval = move.getEvalBranch(traln);
  auto p = toEval.findNodePtr(traln); 
  move.disorientAtNode(traln,p); 
  evaluator.evaluate(traln,toEval, false); 
}


void ExtendedSPR::resetState(TreeAln &traln, PriorBelief &prior )
{
  move.revertTree(traln,prior); 
  // debug_checkTreeConsistency(traln);
  // debug_printTree(traln); 
}


AbstractProposal* ExtendedSPR::clone() const
{
  return new ExtendedSPR( *this); 
}
