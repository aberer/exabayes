#include "BranchLengthMultiplier.hpp"
#include "TreeAln.hpp"
#include "BoundsChecker.hpp"
#include "tune.h"
#include "TreeRandomizer.hpp"
#include "GibbsProposal.hpp"


BranchLengthMultiplier::BranchLengthMultiplier( double _multiplier)
  :  multiplier(_multiplier)
{
  needsFullTraversal = false; 
  this->name = "blMult"; 
  this->category = Category::BRANCH_LENGTHS; 
  relativeWeight = 20;
}

Branch BranchLengthMultiplier::proposeBranch(const TreeAln &traln, Randomness &rand) const 
{
  if(inSetExecution)
    return preparedBranch;
  else  
    return TreeRandomizer::drawBranchUniform(traln, rand); 
}   


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  Branch b = proposeBranch(traln, rand); 
  
  nodeptr p = b.findNodePtr(traln); 
  assert(primaryParameters.size() == 1); 
  auto param = primaryParameters[0].get(); 

  savedBranch = traln.getBranch(p, param); 

  double oldZ = savedBranch.getLength(param);

  double
    drawnMultiplier = 0 ,
    newZ = oldZ; 

  drawnMultiplier= rand.drawMultiplier( multiplier); 
  assert(drawnMultiplier > 0.); 
  newZ = pow( oldZ, drawnMultiplier);
  
  b.setLength(newZ, param); 

  if(not BoundsChecker::checkBranch(b))
    BoundsChecker::correctBranch(b, param); 
  traln.setBranch(b, param); 

  double realMultiplier = log(b.getLength(param)) / log(oldZ); 
  updateHastings(hastings, realMultiplier, name); 

  prior.updateBranchLengthPrior(traln, oldZ, 
				b.getLength(primaryParameters[0].get()), param) ; 
}


void BranchLengthMultiplier::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln) 
{
  evaluator.evaluate(traln,savedBranch, false); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln) 
{
  assert(primaryParameters.size() == 1)  ; 
  traln.setBranch(savedBranch, getPrimaryParameterView()); 
}


void BranchLengthMultiplier::autotune() 
{
  double ratio = sctr.getRatioInLastInterval(); 
  double newParam = tuneParameter(sctr.getBatch(), ratio , multiplier, false);
  multiplier = newParam; 
  sctr.nextBatch();
}


AbstractProposal* BranchLengthMultiplier::clone() const
{
  return new BranchLengthMultiplier(*this);
}


void BranchLengthMultiplier::readFromCheckpointCore(std::istream &in) 
{
  multiplier = cRead<double>(in);
} 

void BranchLengthMultiplier::writeToCheckpointCore(std::ostream &out) const
{
  cWrite(out, multiplier); 
} 


std::pair<Branch,Branch> BranchLengthMultiplier::prepareForSetExecution(TreeAln &traln, Randomness &rand) 
{
  assert(inSetExecution); 
  return std::pair<Branch,Branch>(TreeRandomizer::drawBranchUniform(traln, rand), Branch(0,0)); 
}

