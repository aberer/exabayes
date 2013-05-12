#include "BranchLengthMultiplier.hpp"
#include "eval.h"
#include "TreeAln.hpp"
// #include "proposals.h"
#include "tune.h"

// #define PRINT_MULT


BranchLengthMultiplier::BranchLengthMultiplier(double relativeWeight, double _multiplier)
  :  multiplier(_multiplier)
{
  this->relativeProbability = relativeWeight; 
  this->name = "blMult"; 
  this->category = BRANCH_LENGTHS; 
  ptype = BRANCH_LENGTHS_MULTIPLIER; 
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

  prior.updateBranchLength(branchLengthToReal(tr,oldZ), branchLengthToReal(tr,newZ));

#ifdef PRINT_MULT  
  cout  << setprecision(6) << "bl-multi: " << branchLengthToReal(tr,oldZ) <<   " * " << drawnMultiplier << " = "  << branchLengthToReal(tr, newZ) << endl; 
#endif

  /* according to lakner2008  */
  updateHastings(hastings, drawnMultiplier, name);
 
  /* just doing it for one right here */
  traln.setBranchLengthSave(newZ, 0, p); 
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
  prior.updateBranchLength(newRealZ, oldRealZ); 
  traln.setBranchLengthSave(savedBranch.length[0], 0,p); 
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
