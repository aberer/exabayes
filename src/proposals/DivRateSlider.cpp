#include "DivRateSlider.hpp"

void DivRateSlider::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
}
 
void DivRateSlider::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
}
 
void DivRateSlider::resetState(TreeAln &traln)  
{
}
 
void DivRateSlider::autotune()   
{
}

AbstractProposal* DivRateSlider::clone() const 
{
  return new DivRateSlider(*this); 
}
  
BranchPlain DivRateSlider::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const 
{
}
 
std::vector<nat> DivRateSlider::getInvalidatedNodes(const TreeAln &traln) const 
{
}
  
std::pair<BranchPlain,BranchPlain> DivRateSlider::prepareForSetExecution(TreeAln &traln, Randomness &rand)  
{
}
 
void DivRateSlider::writeToCheckpointCore(std::ostream &out)const   
{
}
  
void DivRateSlider::readFromCheckpointCore(std::istream &in) 
{
}
 
