#include "BranchLengthMultiplier.hpp"



BranchLengthMultiplier::BranchLengthMultiplier(Chain *_chain, double relativeWeight, double _multiplier)
  : chain(_chain), multiplier(_multiplier)
{
  this->relativeProbability = relativeWeight; 
  this->name = "blMult"; 
  this->category = BRANCH_LENGTHS; 
  ptype = BRANCH_LENGTHS_MULTIPLIER; 
}


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  
} 



void BranchLengthMultiplier::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior) 
{
}

void BranchLengthMultiplier::autotune() 
{
  
}


AbstractProposal* BranchLengthMultiplier::clone() const
{
  return new BranchLengthMultiplier(chain, relativeProbability, multiplier);
}
