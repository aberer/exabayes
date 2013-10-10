#include "AminoModelJump.hpp"
#include "Category.hpp"




AminoModelJump::AminoModelJump( vector<aaMatrix_t> matrices)
  : AbstractProposal( Category::AA_MODEL, "aaMat")
{
  relativeWeight = 0.0; // ??? TODO 
}




void AminoModelJump::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator &eval)
{
  assert(NOT_IMPLEMENTED); 
}


void AminoModelJump::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  assert(NOT_IMPLEMENTED); 
}



void AminoModelJump::resetState(TreeAln &traln)
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

