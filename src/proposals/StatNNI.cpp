#include "StatNNI.hpp"
#include "eval.h"
#include "Path.hpp"
#include "output.h"


StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "stNNI" ; 
  this->category = TOPOLOGY; 
  relativeWeight = 5; 
}


void StatNNI::treatOneBranch(nodeptr p, TreeAln &traln, double &hastings, PriorBelief &prior, Randomness &rand)
{
  double m = rand.drawMultiplier( multiplier);
  double oldV = traln.getBranchLength( p,0); 
  double newV = pow(oldV, m); 

  traln.setBranchLengthBounded(newV, 0, p); 
  
  double realM = log(newV) /  log(oldV);   
  updateHastings(hastings, realM, name); 

  auto brPr = secVar[0].getPrior();
  prior.updateBranchLengthPrior(traln,oldV, newV, brPr);
}


void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{  
  tree *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 
  branch b = rand.drawInnerBranchUniform(traln); 

  nodeptr p = findNodeFromBranch(tr, b); 
  
  branch toBePruned = constructBranch(p->next->back->number, p->number ); 
  
  branch switchingBranch; 
  double bl = 0; 
  if(rand.drawRandDouble01() < 0.5)
    {
      switchingBranch = constructBranch(p->back->next->back->number, p->back->number); 
      bl = traln.getBranchLength( p->back->next->back,0); 
    }
  else 
    {
      switchingBranch = constructBranch(p->back->next->next->back->number, p->back->number); 
      bl = traln.getBranchLength( p->back->next->next->back,0); 
    }

  toBePruned.length[0] = traln.getBranchLength( p->next->back,0); 
  b.length[0] = traln.getBranchLength( p,0); 
  switchingBranch.length[0] = bl; 

  path.clear(); 
  path.append( toBePruned);
  path.append( b); 
  path.append( switchingBranch);

  /* switch the branches  */
  nodeptr q = findNodeFromBranch(tr, toBePruned),
    qBack = q->back; 
    
  nodeptr r = findNodeFromBranch(tr, switchingBranch),
    rBack = r->back; 
    
    
#ifdef NNI_MULTIPLY_BL
  bool modifiesBl = false ; 
  for(auto v : secVar)
    modifiesBl |= v.getCategory() == BRANCH_LENGTHS; 

  assert(traln.getNumBranches() == 1 ); 
  if(modifiesBl)
    {
      treatOneBranch(p, traln,hastings, prior, rand);
      treatOneBranch(q, traln,hastings, prior, rand);
      treatOneBranch(r, traln,hastings, prior, rand);
    }
#endif

  traln.clipNode(p,p->back, p->z[0]); 
  traln.clipNode(r, qBack, r->z[0]);
  traln.clipNode(q, rBack, q->z[0]);    
  
  debug_checkTreeConsistency(traln.getTr());
  /* TODO maybe multiply as well */
  
}




void StatNNI::evaluateProposal(TreeAln &traln, PriorBelief &prior)
{
  tree *tr = traln.getTr(); 

  branch b1 = path.at(0), // tobepruned 
    b2  = path.at(2);	// switchingbranch 
  
  branch exchangeBranch = constructBranch(b1.thatNode, b2.thatNode);
  
  nodeptr p = findNodeFromBranch(tr, exchangeBranch),
    q = findNodeFromBranch(tr, invertBranch(exchangeBranch)); 

  newViewGenericWrapper(traln, p, FALSE); 
  newViewGenericWrapper(traln, q, FALSE); 

  evaluateGenericWrapper(traln,p,FALSE ); 
}



void StatNNI::resetState(TreeAln &traln, PriorBelief &prior)  
{
  tree
    *tr = traln.getTr(); 
  int numBranches = traln.getNumBranches(); 

  branch b = path.at(0),
    chosenBranch =path.at(1),
    a = path.at(2); 

  swap(a.thatNode, b.thatNode); 


  /* reset the branches */
  {
    nodeptr q = findNodeFromBranch(tr, a),
      qBack = q->back; 
    nodeptr r = findNodeFromBranch(tr, b), 
      rBack = r->back; 
    
    nodeptr between = findNodeFromBranch(tr, chosenBranch ); 
    
    assert(numBranches == 1); 

    traln.clipNode(between, between->back, chosenBranch.length[0]); 
    traln.clipNode(r, qBack, b.length[0]);  /* r->z */
    traln.clipNode(q, rBack, a.length[0]) ; /* q->z */
  }

  debug_checkTreeConsistency(traln.getTr());

  path.clear();
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
