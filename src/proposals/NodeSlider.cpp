#include "NodeSlider.hpp"
#include "eval/LikelihoodEvaluator.hpp"
#include "system/BoundsChecker.hpp"
#include "priors/AbstractPrior.hpp"

NodeSlider::NodeSlider( double _multiplier)
  : AbstractProposal( Category::BRANCH_LENGTHS, "nodeSlider", 5.,  0,0, false)
  , multiplier(_multiplier)
{
}


std::pair<BranchPlain,BranchPlain> NodeSlider::prepareForSetExecution(TreeAln& traln, Randomness &rand) 
{
  assert(_inSetExecution); 

  auto a = determinePrimeBranch(traln, rand); 
  auto descendents = traln.getDescendents(a); 
  auto b = rand.drawRandDouble01() < 0.5 ? descendents.first : descendents.second; 	
  return std::pair<BranchPlain,BranchPlain>(a,b); 
}


BranchPlain NodeSlider::determinePrimeBranch(const TreeAln &traln, Randomness& rand) const
{
  return TreeRandomizer::drawInnerBranchUniform(traln,rand); 
} 


BranchPlain NodeSlider::proposeBranch(const TreeAln &traln, Randomness &rand) const 
{
  if(_inSetExecution)
    {
      return _preparedBranch;
    }
  else  
    {
      return determinePrimeBranch(traln, rand); 
    }
}


BranchPlain NodeSlider::proposeOtherBranch(const BranchPlain &firstBranch, const TreeAln& traln, Randomness& rand) const 
{
  if(_inSetExecution)
    return _preparedOtherBranch; 
  else 
    {
      auto descendents = traln.getDescendents(firstBranch); 
      return rand.drawRandDouble01() < 0.5 ? descendents.first : descendents.second; 	
    }
}


void NodeSlider::prepareForSetEvaluation( TreeAln &traln, LikelihoodEvaluator& eval) const 
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  
  auto b = BranchPlain(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    

  eval.markDirty( traln, p->number); 
}


void NodeSlider::applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto blParams = getPrimaryParameterView(); 
  auto param = blParams[0];
  assert(blParams.size() == 1); 

  // TODO outer branches have a lower probability of getting drawn? =/ 
  oneBranch = traln.getBranch(proposeBranch(traln, rand),param); 
  otherBranch = traln.getBranch(proposeOtherBranch(oneBranch.toPlain(), traln, rand),param); 

  double oldA = oneBranch.getLength(),
    oldB = otherBranch.getLength();

  double bothZ = oldA * oldB; 

  // correct the multiplier 
  double drawnMultiplier = rand.drawMultiplier(multiplier); 
  double newZ = pow(bothZ, drawnMultiplier); 
  auto testBranch = oneBranch; 
  testBranch.setLength(sqrt(newZ)); 
  if(not BoundsChecker::checkBranch(testBranch))
    {
      BoundsChecker::correctBranch(testBranch); 
      newZ = pow(testBranch.getLength( ),2); 
      drawnMultiplier = log(newZ)   / log(bothZ); 
    }

  double uniScaler = rand.drawRandDouble01(), 
    newA = pow(newZ, uniScaler ),
    newB = pow(newZ, 1-uniScaler); 

  bool problemWithA = false; 
  testBranch.setLength(newA); 
  if(not BoundsChecker::checkBranch(testBranch))
    {
      BoundsChecker::correctBranch(testBranch); 
      newA = testBranch.getLength(); 
      problemWithA = true; 
    }
  testBranch.setLength(newB); 
  if(not BoundsChecker::checkBranch(testBranch))
    {
      BoundsChecker::correctBranch(testBranch); 
      newB = testBranch.getLength(); 
      assert(not problemWithA); // should have been detected, when we checked the multiplier 
    }

  testBranch = oneBranch; 
  testBranch.setLength(newA); 
  // tout << "changing " << oneBranch << " to " << testBranch << std::endl; 
  traln.setBranch(testBranch, param); 
  auto lnPrA = param->getPrior()->getLogProb( ParameterContent{{ testBranch.getInterpretedLength(traln,param) } } )
    /   param->getPrior()->getLogProb( ParameterContent{{ oneBranch.getInterpretedLength(traln,param) } } ); 

  testBranch = otherBranch; 
  testBranch.setLength(newB); 
  // tout << "changing " << otherBranch << " to " << testBranch << std::endl; 
  traln.setBranch(testBranch, param); 
  auto lnPrB = param->getPrior()->getLogProb( ParameterContent{{ testBranch.getInterpretedLength(traln,param) } } )
    /   param->getPrior()->getLogProb( ParameterContent{{ otherBranch.getInterpretedLength(traln,param) } } ); 

  // AbstractProposal::updateHastingsLog(hastings, log(), _name); 
  hastings *= log_double::fromAbs(pow(drawnMultiplier,2));

  prior.addToRatio(lnPrA * lnPrB); 
}

void NodeSlider::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  assert(_primaryParameters.size() == 1); 
  auto b = BranchPlain(middleNode, otherNode); 

  auto parts = _primaryParameters[0]->getPartitions(); 

  for(auto node : getInvalidatedNodes(traln) ) 
    { 
      for(auto part : parts)
	evaluator.markDirty(traln, part, node); 
    }

#ifdef PRINT_EVAL_CHOICE
  tout << "EVAL "  << b << std::endl; 
#endif

  evaluator.evaluatePartitionsWithRoot(traln,b, parts, false); 
}


void NodeSlider::resetState(TreeAln &traln) 
{
  auto view = getPrimaryParameterView();
  assert(view.size() == 1);

  traln.setBranch(oneBranch, view[0]); 
  traln.setBranch(otherBranch, view[0]); 
}


AbstractProposal* NodeSlider::clone() const
{
  return new NodeSlider(*this); 
}  


std::vector<nat> NodeSlider::getInvalidatedNodes(const TreeAln& traln) const
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  return { middleNode } ; 
}
