#include "BranchLengthMultiplier.hpp"
#include "TreeAln.hpp"
#include "tune.h"



BranchLengthMultiplier::BranchLengthMultiplier( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "blMult"; 
  this->category = BRANCH_LENGTHS; 
  relativeWeight = 20;
}


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  tree *tr = traln.getTr(); 
  Branch b =  traln.drawBranchUniform(rand); 
  
  nodeptr p = b.findNodePtr(traln); 
  savedBranch = b; 

  double
    drawnMultiplier = rand.drawMultiplier( multiplier); 
  assert(drawnMultiplier > 0.); 

  assert(traln.getNumBranches() == 1); 

  double oldZ = traln.getBranchLength( p,0);
  savedBranch.setLength( oldZ); 

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


void BranchLengthMultiplier::evaluateProposal(LikelihoodEvaluatorPtr &evaluator,TreeAln &traln, PriorBelief &prior) 
{
  evaluator->evaluate(traln,savedBranch, false); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior) 
{
  nodeptr p = savedBranch.findNodePtr(traln); 
  double tmp = savedBranch.getLength(); 
  traln.setBranchLengthBounded(tmp , 0,p);  
}


void BranchLengthMultiplier::autotune() 
{
  double ratio = sctr.getRatioInLastInterval(); 
  double newParam = tuneParameter(sctr.getBatch(), ratio , multiplier, false);

  // bool up = multiplier  < newParam; 

  // cout << "tuned " << ( up ? " UP ":  " DOWN " )  <<  multiplier << " => " << newParam << "\t ratio=" << setprecision(3) << ratio << endl; 
  
  multiplier = newParam; 
  sctr.nextBatch();
}


AbstractProposal* BranchLengthMultiplier::clone() const
{
  // tout << "cloning "  << name << endl;
  return new BranchLengthMultiplier(*this);
// multiplier
}
