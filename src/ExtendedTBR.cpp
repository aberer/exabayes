#include "ExtendedTBR.hpp"

ExtendedTBR::ExtendedTBR(double _relativeProbability, double extensionProb)
  :  extensionProbability(extensionProb)
{
  name = "eTBR"; 
  relativeProbability = _relativeProbability; 
  category = TOPOLOGY; 
  ptype = E_TBR; 
}



void ExtendedTBR::autotune()
{
  // cannot tune that 
}


void ExtendedTBR::applyToState(TreeAln& traln, PriorManager& prior, double &hastings, Randomness &rand)
{
  
}


void ExtendedTBR::evaluateProposal(TreeAln& traln, PriorManager& prior)
{

}


void ExtendedTBR::resetState(TreeAln &traln, PriorManager& prior)
{
  
}

