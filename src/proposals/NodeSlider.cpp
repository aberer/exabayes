#include "NodeSlider.hpp"
#include "LikelihoodEvaluator.hpp"
#include "BoundsChecker.hpp"

NodeSlider::NodeSlider( double _multiplier)
  : multiplier(_multiplier)
{
  name = "nodeSlider"; 
  this->category = Category::BRANCH_LENGTHS; 
  relativeWeight = 5.; 

  needsFullTraversal = false; 
}


std::pair<BranchPlain,BranchPlain> NodeSlider::prepareForSetExecution(TreeAln& traln, Randomness &rand) 
{
  assert(inSetExecution); 
#ifdef EFFICIENT
  // this is really not okay, we need a generic
  // "prepareforevaluation"method and do not want to care about
  // returned branches
  assert(0); 
#endif
  
  auto a = TreeRandomizer::drawInnerBranchUniform(traln,rand); 
  auto descendents = traln.getDescendents(a); 
  auto b = rand.drawRandDouble01() < 0.5 ? descendents.first : descendents.second; 	
  return std::pair<BranchPlain,BranchPlain>(a,b); 
}


BranchPlain NodeSlider::proposeBranch(const TreeAln &traln, Randomness &rand) const 
{
  if(inSetExecution)
    {
      return preparedBranch;
    }
  else  
    {
      return TreeRandomizer::drawInnerBranchUniform(traln,rand); 
    }
}


BranchPlain NodeSlider::proposeOtherBranch(const BranchPlain &firstBranch, const TreeAln& traln, Randomness& rand) const 
{
  if(inSetExecution)
    return preparedOtherBranch; 
  else 
    {
      auto descendents = traln.getDescendents(firstBranch); 
      return rand.drawRandDouble01() < 0.5 ? descendents.first : descendents.second; 	
    }
}


void NodeSlider::prepareForEvaluation( TreeAln &traln) const 
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  
  auto b = BranchPlain(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 
}


void NodeSlider::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
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
  traln.setBranch(testBranch, param); 
  
  testBranch = otherBranch; 
  testBranch.setLength(newB); 
  traln.setBranch(testBranch, param); 

  updateHastings(hastings, pow(drawnMultiplier,2), name); 

  prior.updateBranchLengthPrior(traln, oldA, newA, param);
  prior.updateBranchLengthPrior(traln, oldB, newB, param);
}

void NodeSlider::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) 
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  
  auto b = BranchPlain(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 

  evaluator.evaluate(traln,b, false); 
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
