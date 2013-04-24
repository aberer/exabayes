#include "nodeSlider.h"


static void insertBranchLength(TreeAln *traln, branch &b)
{
  nodeptr p = findNodeFromBranch(traln->getTr(), b); 
  b.length[0] = traln->getBranchLength( p,0); 
  // getBranchLength( traln->assert(p->number,0) == traln->getBranchLength( p->back->number,0)); 
}

void applyNodeSlider(Chain *chain, proposalFunction *pf)
{
  TreeAln *traln = chain->traln; 
  tree *tr = chain->traln->getTr(); 
  branch oneBranch = drawInnerBranchUniform(chain); 
  insertBranchLength(chain->traln, oneBranch); 
  int numBranch = chain->traln->getNumBranches(); 
  assert(numBranch == 1 ); 

  branch otherBranch;   
  nodeptr p = NULL ; 
  if(isTip(oneBranch.thisNode,tr->mxtips))
    p  = findNodeFromBranch(tr, invertBranch(oneBranch)); 
  else if(isTip(oneBranch.thatNode, tr->mxtips))
    p = findNodeFromBranch(tr, oneBranch); 
  else  
    p = drawRandDouble01(chain) < 0.5  ? findNodeFromBranch(tr, oneBranch) : findNodeFromBranch(tr,invertBranch(oneBranch)); 

  otherBranch.thisNode = p->number; 
  otherBranch.thatNode = 
    drawRandDouble01(chain) < 0.5 ? p->next->back->number : p->next->next->back->number;       

  insertBranchLength(traln, otherBranch); 
  
  // cout << "chose branches " << oneBranch << " and " << otherBranch << endl; 

  stack *thePath = pf->remembrance.modifiedPath; 
  clearStack(thePath);
  pushStack(thePath, oneBranch);
  pushStack(thePath, otherBranch);


  nodeptr nodeA = findNodeFromBranch(tr, oneBranch),
    nodeB = findNodeFromBranch(tr,otherBranch); 

  double bothZ = traln->getBranchLength( nodeA,0) * traln->getBranchLength( nodeB,0); 
  double multiplier = drawMultiplier(chain, pf->parameters.multiplier);
  chain->hastings *= multiplier; 
  double newZ = branchLengthToReal(tr,pow(bothZ,multiplier));  
  double realOldZ = branchLengthToReal(tr, bothZ); 
  chain->hastings *= ( realOldZ / newZ); 
  
  double uniScaler = drawRandDouble01(chain); 
  double aZ = branchLengthToInternal(tr, uniScaler * newZ),
    bZ = branchLengthToInternal(tr, (1-uniScaler) * newZ); 
  
  hookup(nodeA, nodeA->back, &aZ, numBranch); 
  hookup(nodeB, nodeB->back, &bZ, numBranch);  
}


void evaluateNodeSlider(Chain *chain, proposalFunction *pf)
{
  tree *tr = chain->traln->getTr();
  stack *pth = pf->remembrance.modifiedPath; 
  branch otherBranch = getThirdBranch(tr, pth->content[0], pth->content[1]);

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
  newViewGenericWrapper(chain, p, FALSE);  
  evaluateGenericWrapper(chain, p, FALSE );
}



void resetNodeSlider(Chain *chain, proposalFunction *pf)
{
  tree *tr = chain->traln->getTr();
  stack *pth = pf->remembrance.modifiedPath; 
  int numBranches = chain->traln->getNumBranches();

  branch a = popStack(pth),
    b = popStack(pth); 

  nodeptr p = findNodeFromBranch(tr, a),
    q = findNodeFromBranch(tr, b); 
  
  double aZ = a.length[0], 
    bZ = b.length[0]; 
  hookup(p, p->back, &aZ, numBranches); 
  hookup(q,q->back, &bZ, numBranches);   
}


