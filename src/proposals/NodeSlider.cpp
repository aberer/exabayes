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


std::pair<Branch,Branch> NodeSlider::prepareForSetExecution(TreeAln& traln, Randomness &rand) 
{
  assert(inSetExecution); 
#ifdef EFFICIENT
  // this is really not okay, we need a generic
  // "prepareforevaluation"method and do not want to care about
  // returned branches
  assert(0); 
#endif
  
  Branch a = TreeRandomizer::drawInnerBranchUniform(traln,rand); 
  auto descendents = traln.getDescendents(a); 
  Branch b = rand.drawRandDouble01() < 0.5 ? descendents.first : descendents.second; 	
  return std::pair<Branch,Branch>(a,b); 
}


Branch NodeSlider::proposeBranch(const TreeAln &traln, Randomness &rand) const 
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


Branch NodeSlider::proposeOtherBranch(const Branch &firstBranch, const TreeAln& traln, Randomness& rand) const 
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
  
  Branch b(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 
}


void NodeSlider::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  auto blParams = getPrimaryParameterView(); 
  auto param = blParams[0];
  assert(blParams.size() == 1); 

  // TODO outer branches have a lower probability of getting drawn? =/ 
  oneBranch = proposeBranch(traln, rand);
  otherBranch = proposeOtherBranch(oneBranch, traln, rand);

  oneBranch = traln.getBranch(oneBranch,param); 
  otherBranch = traln.getBranch(otherBranch,param); 

  double oldA = oneBranch.getLength(param),
    oldB = otherBranch.getLength(param);

  double bothZ = oldA * oldB; 

  // correct the multiplier 
  double drawnMultiplier = rand.drawMultiplier(multiplier); 
  double newZ = pow(bothZ, drawnMultiplier); 
  auto testBranch = oneBranch; 
  testBranch.setLength(sqrt(newZ), param); 
  if(not BoundsChecker::checkBranch(testBranch))
    {
      BoundsChecker::correctBranch(testBranch, param); 
      newZ = pow(testBranch.getLength( param),2); 
      drawnMultiplier = log(newZ)   / log(bothZ); 
    }

  double uniScaler = rand.drawRandDouble01(), 
    newA = pow(newZ, uniScaler ),
    newB = pow(newZ, 1-uniScaler); 

  bool problemWithA = false; 
  testBranch.setLength(newA, param); 
  if(not BoundsChecker::checkBranch(testBranch))
    {
      BoundsChecker::correctBranch(testBranch, param); 
      newA = testBranch.getLength(param); 
      problemWithA = true; 
    }
  testBranch.setLength(newB, param); 
  if(not BoundsChecker::checkBranch(testBranch))
    {
      BoundsChecker::correctBranch(testBranch, param); 
      newB = testBranch.getLength(param); 
      assert(not problemWithA); // should have been detected, when we checked the multiplier 
    }

  testBranch = oneBranch; 
  testBranch.setLength(newA, param); 
  traln.setBranch(testBranch, param); 
  
  testBranch = otherBranch; 
  testBranch.setLength(newB, param); 
  traln.setBranch(testBranch, param); 

  updateHastings(hastings, pow(drawnMultiplier,2), name); 

  prior.updateBranchLengthPrior(traln, oldA, newA, param);
  prior.updateBranchLengthPrior(traln, oldB, newB, param);
}

void NodeSlider::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) 
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  
  Branch b(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 

  evaluator.evaluate(traln,b, false); 
}


void NodeSlider::resetState(TreeAln &traln) 
{
  auto view = getPrimaryParameterView();
  assert(view.size() == 1); 

  traln.setBranch(oneBranch, view); 
  traln.setBranch(otherBranch, view); 
}


AbstractProposal* NodeSlider::clone() const
{
  return new NodeSlider(*this); 
}  
