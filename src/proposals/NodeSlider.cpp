#include "NodeSlider.hpp"
#include "eval.h"


double NodeSlider::relativeWeight = 5.;


NodeSlider::NodeSlider( double _multiplier)
  : multiplier(_multiplier)
{
  name = "nodeSlider"; 
  this->category = BRANCH_LENGTHS; 
}


static void insertBranchLength(TreeAln &traln, branch &b)
{
  nodeptr p = findNodeFromBranch(traln.getTr(), b); 
  b.length[0] = traln.getBranchLength( p,0); 
}


void NodeSlider::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  tree *tr = traln.getTr(); 
  branch oneBranch = rand.drawInnerBranchUniform(traln); 
  insertBranchLength(traln, oneBranch); 
  int numBranch = traln.getNumBranches(); 
  assert(numBranch == 1 ); 

  branch otherBranch;   
  nodeptr p = NULL ; 
  if(isTip(oneBranch.thisNode,tr->mxtips))
    p  = findNodeFromBranch(tr, invertBranch(oneBranch)); 
  else if(isTip(oneBranch.thatNode, tr->mxtips))
    p = findNodeFromBranch(tr, oneBranch); 
  else  
    p = rand.drawRandDouble01() < 0.5  ? findNodeFromBranch(tr, oneBranch) : findNodeFromBranch(tr,invertBranch(oneBranch)); 

  otherBranch.thisNode = p->number; 
  otherBranch.thatNode = 
    rand.drawRandDouble01() < 0.5 ? p->next->back->number : p->next->next->back->number;       

  insertBranchLength(traln, otherBranch); 

  path.clear();
  path.append(oneBranch);
  path.append(otherBranch);


  nodeptr nodeA = findNodeFromBranch(tr, oneBranch),
    nodeB = findNodeFromBranch(tr,otherBranch); 
  double oldA = traln.getBranchLength( nodeA,0),
    oldB = traln.getBranchLength( nodeB,0);  

  double bothZ = oldA * oldB; 
  double oldBoth = bothZ; 

  double drawnMultiplier = 0,
    newZ = 0;  
  do 
    { 
      drawnMultiplier  =rand.drawMultiplier( multiplier);  
      newZ = pow(bothZ, drawnMultiplier);       
    } while(newZ  < 2 * TreeAln::zMin || newZ > 2 * TreeAln::zMax); 

  double uniScaler = 0, 
    newA = 0,
    newB = 0; 
  do 
    {
      uniScaler = rand.drawRandDouble01();
      newA = pow(bothZ, uniScaler ) ; 
      newB = pow(bothZ, 1-uniScaler); 
    } while (newB < TreeAln::zMin || newA < TreeAln::zMin || newB > TreeAln::zMax || newA > TreeAln::zMax);
  
  traln.setBranchLengthBounded(newA, 0,nodeA);
  traln.setBranchLengthBounded(newB, 0,nodeB); 
  
  updateHastings(hastings, ( log(bothZ) / log(oldBoth)) * drawnMultiplier, name ); 

  auto brPr = randomVariables[0].getPrior(); 

  prior.updateBranchLengthPrior(traln, oldA, newA, brPr);
  prior.updateBranchLengthPrior(traln, oldB, newB, brPr);
}


void NodeSlider::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  tree *tr = traln.getTr();
  branch otherBranch = getThirdBranch(tr, path.at(0), path.at(1));

  nodeptr p = findNodeFromBranch(tr, otherBranch);   
  
  nodeptr q = p->next->back, 
    r = p->next->next->back; 
  
  if(q->x)
    {
      q->x = 0; 
      q->next->x =  1 ; 
    }
  
  if(r->x)
    {
      r->x = 0; 
      r->next->x = 1; 
    }
  // TODO efficient? 
  newViewGenericWrapper(traln, p, FALSE);  
  evaluateGenericWrapper(traln, p, FALSE );
}


void NodeSlider::resetState(TreeAln &traln, PriorBelief &prior) 
{
  tree *tr = traln.getTr();
  int numBranches = traln.getNumBranches();

  assert(path.size() == 2 ); 

  branch &a = path.at(0),
    &b =  path.at(1); 

  nodeptr p = findNodeFromBranch(tr, a),
    q = findNodeFromBranch(tr, b); 
  
  double aZ = a.length[0], 
    bZ = b.length[0]; 
  
  assert(numBranches == 1); 

  traln.clipNode(p, p->back, aZ); 
  traln.clipNode(q,q->back, bZ);   

  path.clear(); 
}


AbstractProposal* NodeSlider::clone() const
{
  return new NodeSlider(*this); 
}  
