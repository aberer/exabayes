#include "axml.h"

#include "TreeLengthMultiplier.hpp"
#include "Randomness.hpp"
#include "TreeAln.hpp"
#include "eval.h"
#include "tune.h"

TreeLengthMultiplier::TreeLengthMultiplier( double _relativeWeight, double _multiplier)
  : multiplier(_multiplier)    
{
  this->relativeProbability = _relativeWeight; 
  this->name = "TL-Mult"; 
  category = BRANCH_LENGTHS; 
  // ptype = TL_MULT; 
}



void TreeLengthMultiplier::multiplyBranchLengthsRecursively(TreeAln& traln, nodeptr p, double multiHere)
{  
  tree *tr = traln.getTr();
  double newZ = pow( traln.getBranchLength( p,0),multiHere); 
  traln.setBranchLengthSave(newZ, 0, p); 

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
#ifdef PRINT_MULT
  cout  << setprecision(6) << "tl-multi with " << rememMultpier << endl; 
#endif

  updateHastings(hastings, rememMultiplier, "TL-Mult");

  initTreeLength =  traln.getTreeLength(); 
  prior.updateBranchLength(initTreeLength, initTreeLength * rememMultiplier);

  multiplyBranchLengthsRecursively(traln , tr->start->back, rememMultiplier); 
}


void TreeLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior)  
{
  tree *tr = traln.getTr();
  multiplyBranchLengthsRecursively(traln, tr->start->back, 1/rememMultiplier);   
  prior.updateBranchLength(initTreeLength * rememMultiplier, initTreeLength);
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
  return new TreeLengthMultiplier(relativeProbability, multiplier);
}  
