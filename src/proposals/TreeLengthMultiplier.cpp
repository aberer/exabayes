#include "axml.h"

#include "TreeLengthMultiplier.hpp"
#include "Randomness.hpp"
#include "TreeAln.hpp"
#include "eval.h"
#include "tune.h"

double TreeLengthMultiplier::relativeWeight = 2 ;

TreeLengthMultiplier::TreeLengthMultiplier( double _multiplier)
  : multiplier(_multiplier)    
{
  this->name = "TL-Mult"; 
  category = BRANCH_LENGTHS; 
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

  auto brPr = primVar[0].getPrior(); 
  
  updateHastings(hastings, rememMultiplier, "TL-Mult");
  initTreeLength = traln.getTreeLengthExpensive(); // TODO? 
  multiplyBranchLengthsRecursively(traln , tr->start->back, rememMultiplier);   
  double newTreeLength = pow( initTreeLength, rememMultiplier ); 
  prior.updateBranchLengthPrior(traln, initTreeLength, newTreeLength, brPr);
}


void TreeLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior)  
{
  tree *tr = traln.getTr();
  multiplyBranchLengthsRecursively(traln, tr->start->back, 1/rememMultiplier);   
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
