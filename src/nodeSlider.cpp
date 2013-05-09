#include "proposalFunction.h"
#include "Randomness.hpp"
#include "Path.hpp"
#include "eval.h"
#include "Chain.hpp"

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
  branch oneBranch = chain->getChainRand()->drawInnerBranchUniform(*traln); 
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
    p = chain->getChainRand()->drawRandDouble01() < 0.5  ? findNodeFromBranch(tr, oneBranch) : findNodeFromBranch(tr,invertBranch(oneBranch)); 

  otherBranch.thisNode = p->number; 
  otherBranch.thatNode = 
    chain->getChainRand()->drawRandDouble01() < 0.5 ? p->next->back->number : p->next->next->back->number;       

  insertBranchLength(traln, otherBranch); 
  
  // cout << "chose branches " << oneBranch << " and " << otherBranch << endl; 

  Path *thePath = pf->remembrance.modifiedPath; 
  // stack *thePath = pf->remembrance.modifiedPath; 
  thePath->clear();
  thePath->append(oneBranch);
  thePath->append(otherBranch);


  nodeptr nodeA = findNodeFromBranch(tr, oneBranch),
    nodeB = findNodeFromBranch(tr,otherBranch); 

  double bothZ = traln->getBranchLength( nodeA,0) * traln->getBranchLength( nodeB,0); 
  double multiplier = chain->getChainRand()->drawMultiplier( pf->parameters.multiplier); 
  chain->hastings *= multiplier; 
  double newZ = branchLengthToReal(tr,pow(bothZ,multiplier));  
  double realOldZ = branchLengthToReal(tr, bothZ); 
#ifdef PRINT_MULT
  cout << setprecision(6) << "nodeslider: " << realOldZ << " * "  << multiplier << " = " << newZ; 
#endif
  chain->hastings *= ( realOldZ / newZ); 
  
  double uniScaler = chain->getChainRand()->drawRandDouble01(); 
  double aZ = branchLengthToInternal(tr, uniScaler * newZ),
    bZ = branchLengthToInternal(tr, (1-uniScaler) * newZ); 
  
  hookup(nodeA, nodeA->back, &aZ, numBranch); 
  hookup(nodeB, nodeB->back, &bZ, numBranch);  
}


void evaluateNodeSlider(Chain *chain, proposalFunction *pf)
{
  tree *tr = chain->traln->getTr();
  Path *pth = pf->remembrance.modifiedPath; 
  branch otherBranch = getThirdBranch(tr, pth->at(0), pth->at(1));

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
  Path *pth = pf->remembrance.modifiedPath; 
  int numBranches = chain->traln->getNumBranches();

  assert(pth->size() == 2 ); 

  branch &a = pth->at(0),
    &b =  pth->at(1); 

  nodeptr p = findNodeFromBranch(tr, a),
    q = findNodeFromBranch(tr, b); 
  
  double aZ = a.length[0], 
    bZ = b.length[0]; 
  hookup(p, p->back, &aZ, numBranches); 
  hookup(q,q->back, &bZ, numBranches);   

  pth->clear(); 
}


