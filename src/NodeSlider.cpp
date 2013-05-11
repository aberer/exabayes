#include "NodeSlider.hpp"
#include "eval.h"


NodeSlider::NodeSlider(double _relativeProbability, double _multiplier)
  : multiplier(_multiplier)
{
  this->relativeProbability = _relativeProbability;   
  name = "nodeSlider"; 
  this->category = BRANCH_LENGTHS; 
  ptype = NODE_SLIDER;   
}


static void insertBranchLength(TreeAln &traln, branch &b)
{
  nodeptr p = findNodeFromBranch(traln.getTr(), b); 
  b.length[0] = traln.getBranchLength( p,0); 
  // getBranchLength( traln.assert(p->number,0) == traln.getBranchLength( p->back->number,0)); 
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

  double bothZ = traln.getBranchLength( nodeA,0) * traln.getBranchLength( nodeB,0); 
  double drawnMultiplier = rand.drawMultiplier( multiplier); 
  hastings *= multiplier; 
  double newZ = branchLengthToReal(tr,pow(bothZ,drawnMultiplier));  
  double realOldZ = branchLengthToReal(tr, bothZ); 
#ifdef PRINT_MULT
  cout << setprecision(6) << "nodeslider: " << realOldZ << " * "  << drawMultiplier << " = " << newZ; 
#endif
  hastings *= ( realOldZ / newZ); 
  
  double uniScaler = rand.drawRandDouble01(); 
  double aZ = branchLengthToInternal(tr, uniScaler * newZ),
    bZ = branchLengthToInternal(tr, (1-uniScaler) * newZ); 
  
  hookup(nodeA, nodeA->back, &aZ, numBranch); 
  hookup(nodeB, nodeB->back, &bZ, numBranch);  
}


void NodeSlider::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{

  tree *tr = traln.getTr();
  // Path *pth = pf->remembrance.modifiedPath; 
  branch otherBranch = getThirdBranch(tr, path.at(0), path.at(1));

  // cout << "third branch is "<< otherBranch << endl; 
  
  nodeptr p = findNodeFromBranch(tr, otherBranch);   
  
  nodeptr q = p->next->back, 
    r = p->next->next->back; 
  
  if(q->x)
    {
      // cout << "disorienting "  << q->number << endl; 
      q->x = 0; 
      q->next->x =  1 ; 
    }
  
  if(r->x)
    {
      // cout << "disorienting "  << r->number << endl; 
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
  hookup(p, p->back, &aZ, numBranches); 
  hookup(q,q->back, &bZ, numBranches);   

  path.clear(); 
}


AbstractProposal* NodeSlider::clone() const
{
  NodeSlider* result = new NodeSlider(relativeProbability, multiplier);
  return result;   
}  
