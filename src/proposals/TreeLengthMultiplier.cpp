#include "axml.h"

#include "TreeLengthMultiplier.hpp"
#include "Randomness.hpp"
#include "TreeAln.hpp"
#include "eval.h"
#include "tune.h"



TreeLengthMultiplier::TreeLengthMultiplier( double _multiplier)
  : multiplier(_multiplier)    
{
  this->name = "TL-Mult"; 
  category = BRANCH_LENGTHS; 
  relativeWeight = 2 ;
}


void TreeLengthMultiplier::multiplyBranchLengthsRecursively(TreeAln& traln, nodeptr p, double multiHere)
{
  tree *tr = traln.getTr();
  double newZ = pow( traln.getBranchLength( p,0),multiHere); 

  traln.setBranchLengthBounded(newZ, 0, p); 

  if(isTip(p->number, tr->mxtips))
    return; 

  nodeptr q = p->next; 
  while(p != q)
    {
      multiplyBranchLengthsRecursively(traln, q->back,multiHere); 
      q = q->next; 
    }
}


void TreeLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  tree *tr = traln.getTr(); 
  rememMultiplier  = rand.drawMultiplier( multiplier);

#ifdef UNSURE
  // for doing this properly, we MUST check, if in any case the
  // multiplier would result in branch lengths that are out of bound
  // and then redraw.

  assert(0); 
#endif

#ifdef EFFICIENT
  // this whole tree length stuff is still highly unsatisfactory 
  assert(0); 
#endif

  auto brPr = primVar[0].getPrior(); 
  
  updateHastings(hastings, pow(rememMultiplier, 2 * traln.getTr()->mxtips - 3  ) , "TL-Mult");

  storedBranches.clear();   
  extractBranches(traln, storedBranches); 
  double initTreeLength = 1,
    newTreeLength = 1;   
  for(auto &b : storedBranches)
    {
      initTreeLength *= b.length[0];
      newTreeLength *= pow(b.length[0], rememMultiplier); 
    }  

  multiplyBranchLengthsRecursively(traln , tr->start->back, rememMultiplier);   
  prior.updateBranchLengthPrior(traln, initTreeLength, newTreeLength, brPr);
}


void TreeLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior)  
{
  for(auto &b : storedBranches)
    {
      nodeptr p = findNodeFromBranch(traln.getTr(), b); 
      traln.clipNode(p,p->back, b.length[0]);
    }
} 


void TreeLengthMultiplier::autotune()
{
  double parameter = multiplier; 

  double newParam = tuneParameter(sctr.getBatch(), sctr.getRatioInLastInterval(), parameter, FALSE);

#ifdef DEBUG_PRINT_TUNE_INFO
  cout << name << ": with ratio " << sctr.getRatioInLastInterval() << ": "<< ((newParam < parameter ) ? "reducing" : "increasing") <<  "\t" << parameter << "," << newParam << endl; 
#endif

  multiplier = newParam; 

  sctr.nextBatch();
}
 
 
void TreeLengthMultiplier::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  evaluateGenericWrapper(traln, traln.getTr()->start, TRUE);
}


AbstractProposal* TreeLengthMultiplier::clone() const
{
  return new TreeLengthMultiplier(*this);
}  
