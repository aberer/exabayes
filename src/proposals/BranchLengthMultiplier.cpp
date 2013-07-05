#include "BranchLengthMultiplier.hpp"
#include "TreeAln.hpp"
#include "BoundsChecker.hpp"
#include "tune.h"
#include "TreeRandomizer.hpp"

#include "GibbsProposal.hpp"


BranchLengthMultiplier::BranchLengthMultiplier( double _multiplier)
  :  multiplier(_multiplier)
{
  this->name = "blMult"; 
  this->category = Category::BRANCH_LENGTHS; 
  relativeWeight = 20;
}

Branch BranchLengthMultiplier::proposeBranch(const TreeAln &traln, Randomness &rand) const 
{
  return TreeRandomizer::drawBranchUniform(traln, rand); 
}   


void BranchLengthMultiplier::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  Branch b = proposeBranch(traln, rand); 
  
  nodeptr p = b.findNodePtr(traln); 
  savedBranch = b; 

  double oldZ = traln.getBranch(p).getLength();
  savedBranch.setLength( oldZ); 
  assert(traln.getNumBranches() == 1); 
  
  double
    drawnMultiplier = 0 ,
    newZ = oldZ; 
  
  do 
    {
      drawnMultiplier= rand.drawMultiplier( multiplier); 
      assert(drawnMultiplier > 0.); 
      newZ = pow( oldZ, drawnMultiplier);
    } while(newZ < BoundsChecker::zMin || BoundsChecker::zMax < newZ ) ; 

  // TODO 
 // double newZ = oldZ; 
 //  if(not traln.isCollapsed(b))    

 //  else 
 //    dra wnMultiplier = 1; 

  traln.setBranch(Branch(p->number, p->back->number, newZ)); 

  double realMultiplier = log(newZ) / log(oldZ); 
  updateHastings(hastings, realMultiplier, name); 

  // cout << "acessing "<< *(primVar[0] )   << endl; 
  auto relPrior =  primVar[0]->getPrior(); 
  prior.updateBranchLengthPrior(traln, oldZ, newZ,relPrior) ; 
}


void BranchLengthMultiplier::evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, PriorBelief &prior) 
{
  evaluator.evaluate(traln,savedBranch, false); 
}

 
void BranchLengthMultiplier::resetState(TreeAln &traln, PriorBelief &prior) 
{
  traln.setBranch(savedBranch); 
}


void BranchLengthMultiplier::autotune() 
{
  double ratio = sctr.getRatioInLastInterval(); 
  double newParam = tuneParameter(sctr.getBatch(), ratio , multiplier, false);

  // bool up = multiplier  < newParam; 

  // cout << "tuned " << ( up ? " UP ":  " DOWN " )  <<  multiplier << " => " << newParam << "\t ratio=" << setprecision(3) << ratio << endl; 
  
  multiplier = newParam; 
  sctr.nextBatch();
}


AbstractProposal* BranchLengthMultiplier::clone() const
{
  // tout << "cloning "  << name << endl;
  return new BranchLengthMultiplier(*this);
// multiplier
}
