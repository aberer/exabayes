#include "BranchLengthMultiplier.hpp"
#include "eval.h"
#include "TreeAln.hpp"
#include "tune.h"

BranchLengthMultiplier::BranchLengthMultiplier(double relativeWeight, double _multiplier)
  :  multiplier(_multiplier)
{
  this->relativeProbability = relativeWeight; 
  this->name = "blMult"; 
  this->category = BRANCH_LENGTHS; 
  // ptype = BRANCH_LENGTHS_MULTIPLIER; 
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
  
  /* TODO how do we do that wiht multiple bls per branch?  */
  assert(traln.getNumBranches() == 1); 

  double oldZ = traln.getBranchLength( p,0);
  savedBranch.length[0] = oldZ; 
  double newZ = pow( oldZ, drawnMultiplier) ; 


#ifdef PRINT_MULT  
  cout  << setprecision(6) << name << branchLengthToReal(tr,oldZ) <<   " * " << drawnMultiplier << " = "  << branchLengthToReal(tr, newZ) << endl; 
#endif

  /* just doing it for one right here */
  traln.setBranchLengthBounded(newZ, 0, p); 

  double realOld = branchLengthToReal(tr,oldZ); 
  double realNew = branchLengthToReal(tr,newZ); 
  prior.updateBranchLength(realOld, realNew);
  drawnMultiplier =  realNew / realOld ; 

  /* according to lakner2008  */
  updateHastings(hastings, drawnMultiplier, name);
} 



void BranchLengthMultiplier::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = findNodeFromBranch(traln.getTr(), savedBranch); 
  evaluateGenericWrapper(traln,p, FALSE ); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = findNodeFromBranch(traln.getTr(), savedBranch); 
  double newRealZ = branchLengthToReal(traln.getTr(), p->z[0]),
    oldRealZ = branchLengthToReal(traln.getTr(), savedBranch.length[0]); 
  traln.setBranchLengthBounded(savedBranch.length[0], 0,p); 
  prior.updateBranchLength(newRealZ, oldRealZ); 
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
  return new BranchLengthMultiplier(relativeProbability, multiplier);
}
