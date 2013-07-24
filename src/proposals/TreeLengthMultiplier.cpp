#include "axml.h"

#include "BoundsChecker.hpp"
#include "TreeLengthMultiplier.hpp"
#include "Randomness.hpp"
#include "TreeAln.hpp"
#include "tune.h"



TreeLengthMultiplier::TreeLengthMultiplier( double _multiplier)
  : multiplier(_multiplier)    
{
  this->name = "TL-Mult"; 
  category = Category::BRANCH_LENGTHS; 
  relativeWeight = 2 ;
}


void TreeLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  storedBranches.clear(); 
  storedBranches = traln.extractBranches();

  std::vector<Branch> newBranches = storedBranches; 
  
  double treeScaler = rand.drawMultiplier(multiplier); 
  double initTL = 1; 
  double newTL = 1; 
  // std::cout << "drew " <<  treeScaler << std::endl; 

  for(auto &b : newBranches)
    {
      auto initLength = b.getLength();
      initTL *= initLength; 

      b.setLength(  pow(initLength, treeScaler) ); 
      
      if( not BoundsChecker::checkBranch(b))
	{
	  // std::cout << "correction!" << std::endl; 
	  BoundsChecker::correctBranch(b);
	}
      
      double realScaling = log(b.getLength())  /  log(initLength); 
      updateHastings(hastings, realScaling, "TL-multi"); 
      newTL *= b.getLength(); 
    }

  for(auto &b : newBranches)
    traln.setBranch(b);

  auto brPr = primVar[0]->getPrior(); 
  prior.updateBranchLengthPrior(traln, initTL, newTL, brPr);
}


void TreeLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior)  
{
  for(auto &b : storedBranches)
    traln.setBranch(b); 
} 


void TreeLengthMultiplier::autotune()
{
  double parameter = multiplier; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, FALSE);

#ifdef DEBUG_PRINT_TUNE_INFO
  std::cout << name << ": with ratio " << sctr.getRatioInLastInterval()
	    << ": "<< ((newParam < parameter ) ? "reducing" : "increasing") <<  "\t" << parameter << "=>" << newParam << std::endl; 
#endif

  multiplier = newParam; 

  sctr.nextBatch();
}
 
 
void TreeLengthMultiplier::evaluateProposal(  LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior) 
{
  evaluator->evaluate(traln,Branch(traln.getTr()->start->number,traln.getTr()->start->back->number), true); 
}


AbstractProposal* TreeLengthMultiplier::clone() const
{
  return new TreeLengthMultiplier(*this);
}  


void TreeLengthMultiplier::readFromCheckpointCore(std::ifstream &in)
{
  multiplier = cRead<double>(in); 
} 

void TreeLengthMultiplier::writeToCheckpointCore(std::ofstream &out) 
{
  cWrite(out, multiplier); 
} 


