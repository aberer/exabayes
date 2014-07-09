#include "DivTimeProposal.hpp"


DivTimeProposal::DivTimeProposal(   )  
  : AbstractProposal(Category::BRANCH_LENGTHS, "divTimeDirich", 10., 0,0, true)
    // TODO: should be a 
{
  // todo initialize parameters, e.g. a parameter to tune the dirichlet proposal
} 


void DivTimeProposal::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
  assert(0); 
} 

void DivTimeProposal::evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion)
{
  // trivial full evaluaet (just chekc the toher proposals )
  assert(0); 
} 


void DivTimeProposal::resetState(TreeAln &traln) 
{
  // reset the previous state of traln 
  assert(0); 
}



void DivTimeProposal::autotune() 
{
  // tune proposal parameter? 
}


AbstractProposal* DivTimeProposal::clone() const
{
  return new DivTimeProposal(*this); 
}  


BranchPlain DivTimeProposal::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const
{
  // TODO implement, if you want to have proposal sets (not necessary iniitiallyz)
  assert(0); 
} 


std::vector<nat> DivTimeProposal::getInvalidatedNodes(const TreeAln &traln) const
{
  // TODO implement, if you want to have proposal sets (not necessary iniitiallyz)
  assert(0); 
} 



std::pair<BranchPlain,BranchPlain> DivTimeProposal::prepareForSetExecution(TreeAln &traln, Randomness &rand) 
{
  // TODO implement, if you want to have proposal sets (not necessary iniitiallyz)
  assert(0); 
} 


void DivTimeProposal::writeToCheckpointCore(std::ostream &out)const
{
  // todo implemeent, in case you have proposal parameters (e.g., alpha of dirichlet)
  // see other proposals 
}


void DivTimeProposal::readFromCheckpointCore(std::istream &in)
{
  // todo implemeent, in case you have proposal parameters (e.g., alpha of dirichlet)
  // see other proposals 
} 
