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
  needsFullTraversal = true; 
}


// TreeLengthMultiplier::TreeLengthMultiplier(const TreeLengthMultiplier& rhs)
//   : AbstractProposal(rhs)
// {
// }


void TreeLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  storedBranches.clear(); 
  
  auto blParam = primaryParameters[0].get(); 

  assert(primaryParameters.size() == 1); 

  storedBranches = traln.extractBranches(blParam);

  std::vector<Branch> newBranches = storedBranches; 
  
  double treeScaler = rand.drawMultiplier(multiplier); 
  double initTL = 1; 
  double newTL = 1; 
  // std::cout << "drew " <<  treeScaler << std::endl; 

  for(auto &b : newBranches)
    {
      auto initLength = b.getLength( blParam);
      initTL *= initLength; 

      b.setLength(  pow(initLength, treeScaler), blParam ); 
      
      if( not BoundsChecker::checkBranch(b))
	{
	  // std::cout << "correction!" << std::endl; 
	  BoundsChecker::correctBranch(b, blParam);
	}
      
      double realScaling = log(b.getLength(blParam)) / log(initLength); 
      updateHastings(hastings, realScaling, "TL-multi"); 
      newTL *= b.getLength(blParam); 
    }

  for(auto &b : newBranches)
    {
      std::vector<AbstractParameter*> tmp = {blParam}; 
      traln.setBranch(b, tmp);
    }

  prior.updateBranchLengthPrior(traln, initTL, newTL, blParam);
}


void TreeLengthMultiplier::resetState(TreeAln &traln)  
{
  for(auto &b : storedBranches)
    traln.setBranch(b, getPrimaryParameterView()); 
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
 
 
void TreeLengthMultiplier::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) 
{
  evaluator.evaluate(traln,Branch(traln.getTr()->start->number,traln.getTr()->start->back->number), true); 
}


AbstractProposal* TreeLengthMultiplier::clone() const
{
  return new TreeLengthMultiplier(*this);
}  


void TreeLengthMultiplier::readFromCheckpointCore(std::istream &in)
{
  multiplier = cRead<double>(in); 
} 

void TreeLengthMultiplier::writeToCheckpointCore(std::ostream &out) const
{
  cWrite(out, multiplier); 
} 


