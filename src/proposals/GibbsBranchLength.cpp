#include "GibbsBranchLength.hpp"


GibbsBranchLength::GibbsBranchLength()
  : BranchLengthMultiplier(0)
{
  name = "estGibbsBL"; 
  category = Category::BRANCH_LENGTHS; 
  relativeWeight = 20 ; 
} 



GibbsBranchLength::GibbsBranchLength(const GibbsBranchLength& rhs) 
  : BranchLengthMultiplier(rhs)
{
}


void GibbsBranchLength::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) 
{
  Branch b = proposeBranch(traln, rand);     
  
  auto params = getBranchLengthsParameterView();
  assert(params.size() == 1 );     
  auto param = params[0]; 

  b  = traln.getBranch(b,param); 
  double initBl = b.getLength(param); 
  savedBranch = b; 

  auto newBranch = GibbsProposal::drawFromEsitmatedPosterior(b, eval, traln, rand,  MAX_ITER, hastings, param); 
  traln.setBranch(newBranch, param); 

  prior.updateBranchLengthPrior(traln,initBl, newBranch.getLength(param), param);  
}
