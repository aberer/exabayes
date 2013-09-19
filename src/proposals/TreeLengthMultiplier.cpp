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


void TreeLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  storedBranches.clear(); 
  
  auto blParam = primaryParameters[0].get(); 

  assert(primaryParameters.size() == 1); 

  storedBranches = traln.extractBranches(blParam);

  auto  newBranches = storedBranches; 
  
  double treeScaler = rand.drawMultiplier(multiplier); 
  double initTL = 0; 
  double newTL = 0; 
  // std::cout << "drew " <<  treeScaler << std::endl; 

  for(auto &b : newBranches)
    {
      auto initLength = b.getLength();
      initTL += log(initLength); 

      b.setLength(  pow(initLength, treeScaler) ); 
      
      if( not BoundsChecker::checkBranch(b))
	BoundsChecker::correctBranch(b);

      double realScaling = log(b.getLength()) / log(initLength); 
      updateHastings(hastings, realScaling, "TL-multi"); 
      newTL += log(b.getLength()); 
    }

  for(auto &b : newBranches)
    traln.setBranch(b, blParam);

  // well ... =/ 
  initTL = exp(initTL); 
  newTL = exp(newTL); 

  // tout  << MAX_SCI_PRECISION << "treeScaler=" << treeScaler << "\tinitTL=" << initTL << "\tnewTL="  << newTL << std::endl; 

  prior.updateBranchLengthPrior(traln, initTL, newTL, blParam);
}


void TreeLengthMultiplier::resetState(TreeAln &traln)  
{
  auto blParams = getPrimaryParameterView(); 
  assert(blParams.size( )== 1 ) ;
  auto blParam = blParams[0]; 
  for(auto &b : storedBranches)
    traln.setBranch(b, blParam); 
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
  auto b = BranchPlain(traln.getTr()->start->number,traln.getTr()->start->back->number); 
  evaluator.evaluate(traln,b, true); 
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


