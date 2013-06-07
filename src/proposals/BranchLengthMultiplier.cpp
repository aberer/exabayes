#include "BranchLengthMultiplier.hpp"
#include "eval.h"
#include "TreeAln.hpp"
#include "tune.h"

double BranchLengthMultiplier::relativeWeight = 20;

BranchLengthMultiplier::BranchLengthMultiplier( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "blMult"; 
  this->category = BRANCH_LENGTHS; 
}


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  tree *tr = traln.getTr(); 
  branch b =  rand.drawBranchUniform(traln); 
  
  nodeptr p = findNodeFromBranch(tr, b); 
  savedBranch = b; 

  double
    drawnMultiplier = rand.drawMultiplier( multiplier); 
  assert(drawnMultiplier > 0.); 

  assert(traln.getNumBranches() == 1); 

  double oldZ = traln.getBranchLength( p,0);
  savedBranch.length[0] = oldZ; 

  double newZ = oldZ; 
  if(not traln.isCollapsed(b))    
    newZ = pow( oldZ, drawnMultiplier);
  else 
    drawnMultiplier = 1; 

  /* just doing it for one right here */
  traln.setBranchLengthBounded(newZ, 0, p); 

  double realMultiplier = log(newZ) / log(oldZ); 
  updateHastings(hastings, realMultiplier, name); 

  prior.updateBranchLengthPrior(traln, oldZ, newZ, primVar[0].getPrior()); 
} 


void BranchLengthMultiplier::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = findNodeFromBranch(traln.getTr(), savedBranch); 
  evaluateGenericWrapper(traln,p, FALSE ); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = findNodeFromBranch(traln.getTr(), savedBranch); 
  traln.setBranchLengthBounded(savedBranch.length[0], 0,p);  
}


void BranchLengthMultiplier::autotune() 
{
  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), multiplier, FALSE);
  multiplier = newParam; 
  sctr.nextBatch();
}


AbstractProposal* BranchLengthMultiplier::clone() const
{
  // tout << "cloning "  << name << endl;
  return new BranchLengthMultiplier(*this);
// multiplier
}
