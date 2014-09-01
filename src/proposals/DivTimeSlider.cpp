#include "DivTimeSlider.hpp"

void DivTimeSlider::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
}


void DivTimeSlider::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  assert(0); 
} 


void DivTimeSlider::resetState(TreeAln &traln)  
{
}
 


void DivTimeSlider::autotune()   
{
}

AbstractProposal* DivTimeSlider::clone() const 
{
  return new DivTimeSlider(*this); 
}


BranchPlain DivTimeSlider::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
{

}

 
std::vector<nat> DivTimeSlider::getInvalidatedNodes(const TreeAln &traln) const 
{
  
}


std::pair<BranchPlain,BranchPlain> DivTimeSlider::prepareForSetExecution(TreeAln &traln, Randomness &rand)  
{
}
 
void DivTimeSlider::writeToCheckpointCore(std::ostream &out)const   
{
}
  
void DivTimeSlider::readFromCheckpointCore(std::istream &in) 
{
}
 
