#include "StatNNI.hpp"
#include "eval.h"
#include "Path.hpp"
#include "output.h"


StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "stNNI" ; 
  this->category = TOPOLOGY; 
  relativeWeight = 5; 
}


void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{    
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 

  Branch b( rand.drawInnerBranchUniform(traln) ) ; 
  nodeptr p = b.findNodeFromBranch(traln); 

  pair<int,int> switching(p->next->back->number,
			  rand.drawRandDouble01() < 0.5  
			  ? p->back->next->back->number
			  : p->back->next->next->back->number
			  ); 

  move.init( traln, b,  switching);

  vector<shared_ptr<AbstractPrior> > priors; 
  bool multiplyBranches = false; 
  for(auto v : secVar)
    {
      multiplyBranches |= v.getCategory() == BRANCH_LENGTHS; 
      priors.push_back(v.getPrior()) ; 
#ifdef UNSURE
      // fixed prior? 
      assert(0); 
#endif
    }

#ifdef NO_SEC_BL_MULTI
  multiplyBranches = false; 
#endif

  if(multiplyBranches)
    move.multiplyAllBranches(traln, hastings, prior, rand, multiplier, priors, name);

  move.apply(traln);
  debug_checkTreeConsistency(traln);
}




void StatNNI::evaluateProposal(TreeAln &traln, PriorBelief &prior)
{
  move.disortient(traln); 
  nodeptr p = move.getEvalBranch().findNodeFromBranch(traln);
  evaluateGenericWrapper(traln, p,FALSE); 
}



void StatNNI::resetState(TreeAln &traln, PriorBelief &prior)  
{
  move.revert(traln); 
  debug_checkTreeConsistency(traln);
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
