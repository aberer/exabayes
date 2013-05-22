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

  double bothZ = traln.getBranchLength( nodeA,0) * traln.getBranchLength( nodeB,0); 
  double drawnMultiplier = rand.drawMultiplier( multiplier); 

  double newZ = branchLengthToReal(tr,pow(bothZ,drawnMultiplier));
#ifdef PRINT_MULT
  cout << setprecision(6) << name << realOldZ << " * "  << drawnMultiplier << " = " << newZ << endl;  
#endif

  double uniScaler = rand.drawRandDouble01(); 
  double aZ = branchLengthToInternal(tr, uniScaler * newZ),
    bZ = branchLengthToInternal(tr, (1-uniScaler) * newZ); 

#if TODO 
  double zABefore = branchLengthToReal(tr, traln.getBranchLength(nodeA,0)),
    zBBefore =branchLengthToReal(tr,traln.getBranchLength(nodeB, 0)); 
#endif

  traln.clipNode(nodeA, nodeA->back, aZ); 
  traln.clipNode(nodeB, nodeB->back, bZ); 

  // update the hastings with the REAL value of the branch lengths
  // after multiplication. This has to be so complicated, since
  // clipNode sets the bl in a bounded region.
  updateHastings(hastings,   branchLengthToReal(tr, aZ * bZ ) / branchLengthToReal(tr, bothZ), "nodeSlider" ); 
  // TODO@kassian: is this hastings correct? 

#ifdef UNSURE
  assert(0); 
#endif

#if TODO  
  prior.updateBranchLength( zABefore , branchLengthToReal(traln.getTr(), aZ));
  prior.updateBranchLength( zBBefore , branchLengthToReal(traln.getTr(), bZ)); 
#endif
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
#if TODO 
  prior.updateBranchLength(branchLengthToReal(traln.getTr(), p->z[0]), branchLengthToReal(traln.getTr(), aZ));
  prior.updateBranchLength(branchLengthToReal(traln.getTr(), q->z[0]), branchLengthToReal(traln.getTr(), bZ));
#endif
  traln.clipNode(p, p->back, aZ); 
  traln.clipNode(q,q->back, bZ);   

  path.clear(); 
}


AbstractProposal* NodeSlider::clone() const
{
  return new NodeSlider(*this); 
  // NodeSlider* result = new NodeSlider( multiplier);
  // return result;   
}  
