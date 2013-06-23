#include "NodeSlider.hpp"
#include "eval.h"
#include "LikelihoodEvaluator.hpp"

NodeSlider::NodeSlider( double _multiplier)
  : multiplier(_multiplier)
{
  name = "nodeSlider"; 
  this->category = BRANCH_LENGTHS; 
  relativeWeight = 5.; 
}


void NodeSlider::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  oneBranch = traln.drawInnerBranchUniform(rand);     
  if(rand.drawRandDouble01() < 0.5 )
    oneBranch.invert();

  nodeptr p = oneBranch.findNodePtr(traln); 
  nodeptr q = rand.drawRandDouble01() < 0.5  ? p->next->back : p->next->next->back; 
  oneBranch.setLength(traln.getBranchLength(p,0)); 
  otherBranch = Branch(p->number, q->number, traln.getBranchLength(q,0)); 
  
  double oldA = oneBranch.getLength(),
    oldB = otherBranch.getLength();

  double bothZ = oldA * oldB; 
  double oldBoth = bothZ; 

  double drawnMultiplier = 0,
    newZ = 0;  

  drawnMultiplier  =rand.drawMultiplier( multiplier);  
  newZ = pow(bothZ, drawnMultiplier);       
  
  double uniScaler = 0, 
    newA = 0,
    newB = 0; 
  uniScaler = rand.drawRandDouble01();
  newA = pow(newZ, uniScaler ) ; 
  newB = pow(newZ, 1-uniScaler); 

  traln.clipNode(oneBranch.findNodePtr(traln), 
		oneBranch.getInverted().findNodePtr(traln),
		newA); 
  
  traln.clipNode(otherBranch.findNodePtr(traln),
		 otherBranch.getInverted().findNodePtr(traln),
		 newB); 

#ifdef UNSURE
  assert(0); 
  // extremely unsure about all that 
#endif

  double realBoth = newA * newB; 
  
  double realMulti = realBoth / oldBoth ; 
  
  updateHastings(hastings, ( log(realBoth) / log(oldBoth)) * realMulti, name ); 

  auto brPr = primVar[0].getPrior(); 

  prior.updateBranchLengthPrior(traln, oldA, newA, brPr);
  prior.updateBranchLengthPrior(traln, oldB, newB, brPr);
}

void NodeSlider::evaluateProposal(  LikelihoodEvaluatorPtr &evaluator, TreeAln &traln, PriorBelief &prior) 
{
  nat middleNode = oneBranch.getCommonNode(otherBranch); 
  nat otherNode = oneBranch.getOtherNode(middleNode); 
  
  Branch b(middleNode, otherNode); 
  nodeptr p = b.findNodePtr(traln);    
  LikelihoodEvaluator::disorientNode( p); 

  evaluator->evaluate(traln,b, false); 
}


void NodeSlider::resetState(TreeAln &traln, PriorBelief &prior) 
{
  oneBranch.applyToTree(traln); 
  otherBranch.applyToTree(traln); 
}


AbstractProposal* NodeSlider::clone() const
{
  return new NodeSlider(*this); 
}  
