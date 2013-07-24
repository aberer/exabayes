#include "NodeSlider.hpp"
#include "LikelihoodEvaluator.hpp"
#include "BoundsChecker.hpp"

NodeSlider::NodeSlider( double _multiplier)
  : multiplier(_multiplier)
{
  name = "nodeSlider"; 
  this->category = Category::BRANCH_LENGTHS; 
  relativeWeight = 5.; 
}


void NodeSlider::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  // TODO outer branches have a lower probability of getting drawn? =/ 
  oneBranch = TreeRandomizer::drawInnerBranchUniform(traln,rand); 

  nodeptr p = oneBranch.findNodePtr(traln); 
  // TODO correct here for the problem mentioned above   
  nodeptr q = rand.drawRandDouble01() < 0.5  ? p->next->back : p->next->next->back; 
  oneBranch.setLength(traln.getBranch(p).getLength()); 
  otherBranch = Branch(p->number, q->number, traln.getBranch(q).getLength()); 
  
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
      newZ = pow(testBranch.getLength(),2); 
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

  traln.clipNode(oneBranch.findNodePtr(traln), 
		oneBranch.getInverted().findNodePtr(traln),
		newA); 
  
  traln.clipNode(otherBranch.findNodePtr(traln),
		 otherBranch.getInverted().findNodePtr(traln),
		 newB); 
  updateHastings(hastings, drawnMultiplier * drawnMultiplier, name); 
  auto brPr = primVar[0]->getPrior(); 

  prior.updateBranchLengthPrior(traln, oldA, newA, brPr);
  prior.updateBranchLengthPrior(traln, oldB, newB, brPr);
}

void NodeSlider::evaluateProposal(  LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior) 
{
  nat middleNode = oneBranch.getIntersectingNode(otherBranch) ;
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  
  Branch b(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 

  evaluator->evaluate(traln,b, false); 
}


void NodeSlider::resetState(TreeAln &traln, PriorBelief &prior) 
{
  traln.setBranch(oneBranch); 
  traln.setBranch(otherBranch); 
}


AbstractProposal* NodeSlider::clone() const
{
  return new NodeSlider(*this); 
}  
