#include "WrappedProposal.hpp"


WrappedProposal::WrappedProposal(proposalFunction *pf, Chain *_chain)
  : pfun(pf)
  , chain(_chain)
{
  name = pf->name; 
  category = pf->category; 
  relativeProbability = pf->relativeWeight; 
  ptype = pf->ptype; 
}

void WrappedProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  pfun->apply_func(chain, pfun);
}


void WrappedProposal::evaluateProposal(TreeAln &traln, PriorBelief &prior)
{
  pfun->eval_lnl(chain, pfun);
}

void WrappedProposal::resetState(TreeAln &traln, PriorBelief &prior) 
{
  pfun->reset_func(chain, pfun); 
}


void WrappedProposal::autotune()
{
  if(pfun->autotune != NULL)
    pfun->autotune(pfun, &sctr);
}


WrappedProposal* WrappedProposal::clone() const
{  
  WrappedProposal* result = new WrappedProposal(pfun, chain);
  return result; 
} 
