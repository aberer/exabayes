#include "AminoModelJump.hpp"


AminoModelJump::AminoModelJump(double _relativeProbability, vector<aaMatrix_t> matrices)
{
  this->relativeProbability = _relativeProbability; //  the constructor argument should not have the exact same name as the member variable
  name = "aaMat"; 	
  category = SUBSTITUTION_RATES ; 	// check out categoryType.h

}




void AminoModelJump::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  assert(NOT_IMPLEMENTED); 
}


void AminoModelJump::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  assert(NOT_IMPLEMENTED); 
}



void AminoModelJump::resetState(TreeAln &traln, PriorBelief &prior)
{  
  assert(NOT_IMPLEMENTED); 
}



void AminoModelJump::autotune() 
{
  // do nothing
}


AbstractProposal* AminoModelJump::clone() const 
{
  return new AminoModelJump(relativeProbability, matrices);
}

