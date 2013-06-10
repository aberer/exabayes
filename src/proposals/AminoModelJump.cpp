#include "AminoModelJump.hpp"




AminoModelJump::AminoModelJump(vector<aaMatrix_t> matrices)
{
  name = "aaMat"; 	
  category = AA_MODEL ; 
  relativeWeight = 0.0; // ??? TODO 
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
  // return new AminoModelJump( matrices);
  return new AminoModelJump(*this); 
}

