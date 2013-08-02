#include "StatNNI.hpp"
#include "Path.hpp"
#include "TreeRandomizer.hpp"


StatNNI::StatNNI( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "stNNI" ; 
  this->category = Category::TOPOLOGY; 
  relativeWeight = 5; 
}


void StatNNI::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{    
  int numBranches = traln.getNumBranches();
  assert(numBranches == 1); 

  Branch b( TreeRandomizer::drawInnerBranchUniform(traln, rand) ) ; 
  nodeptr p = b.findNodePtr(traln); 

  Branch switchingBranch = Branch( 
				  rand.drawRandDouble01() < 0.5  
				  ? p->back->next->back->number
				  : p->back->next->next->back->number, 
				  p->back->number
				   ); 
  
  move.extractMoveInfo(traln, {b, switchingBranch}); 

  std::vector<AbstractPrior* > priors; 
  bool multiplyBranches = false; 
  for(auto &v : secondaryParameters)
    {
      multiplyBranches |= v->getCategory() == Category::BRANCH_LENGTHS; 
      priors.push_back(v->getPrior()) ; 
#ifdef UNSURE
      // fixed prior? 
      assert(0); 
#endif
    }

#ifdef NO_SEC_BL_MULTI
  multiplyBranches = false; 
#endif

  if(multiplyBranches)
    move.multiplyBranches(traln, rand, hastings, prior, multiplier, priors); 

  move.applyToTree(traln);
  // debug_checkTreeConsistency(traln);
}




void StatNNI::evaluateProposal(  LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior)
{
  Branch evalBranch = move.getEvalBranch(traln); 
  nodeptr p = evalBranch.findNodePtr(traln);
  move.disorientAtNode(traln, p);
  evaluator->evaluate(traln, evalBranch, false); 
}



void StatNNI::resetState(TreeAln &traln, PriorBelief &prior)  
{
  move.revertTree(traln,prior); 
  // debug_checkTreeConsistency(traln);
}


AbstractProposal* StatNNI::clone()  const
{
  return new StatNNI( *this );
}
