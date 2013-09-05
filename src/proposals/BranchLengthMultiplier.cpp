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

BranchPlain BranchLengthMultiplier::proposeBranch(const TreeAln &traln, Randomness &rand) const 
{
  if(inSetExecution)
    return preparedBranch;
  else  
    return TreeRandomizer::drawBranchUniform(traln, rand); 
}   


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto b = proposeBranch(traln, rand).toBlDummy(); 
  
  // nodeptr p = b.findNodePtr(traln); 

  assert(primaryParameters.size() == 1); 
  auto param = primaryParameters[0].get(); 

  savedBranch = traln.getBranch(b.toPlain(), param); 

  double oldZ = savedBranch.getLength();

  double
    drawnMultiplier = 0 ,
    newZ = oldZ; 

  drawnMultiplier= rand.drawMultiplier( multiplier); 
  assert(drawnMultiplier > 0.); 
  newZ = pow( oldZ, drawnMultiplier);
  
  b.setLength(newZ); 

  if(not BoundsChecker::checkBranch(b))
    BoundsChecker::correctBranch(b); 
  traln.setBranch(b, param); 

  double realMultiplier = log(b.getLength()) / log(oldZ); 
  updateHastings(hastings, realMultiplier, name); 

  prior.updateBranchLengthPrior(traln, oldZ, b.getLength(), param) ; 
}


void BranchLengthMultiplier::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln) 
{
  evaluator.evaluate(traln,savedBranch.toPlain(), false); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln) 
{
  assert(primaryParameters.size() == 1)  ; 
  auto params = getBranchLengthsParameterView(); 
  assert(params.size() == 1); 
  auto param = params[0]; 
  traln.setBranch(savedBranch, param); 
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


std::pair<BranchPlain,BranchPlain> BranchLengthMultiplier::prepareForSetExecution(TreeAln &traln, Randomness &rand) 
{
  assert(inSetExecution); 
  return std::pair<BranchPlain,BranchPlain>(TreeRandomizer::drawBranchUniform(traln, rand), BranchPlain(0,0)); 
}

