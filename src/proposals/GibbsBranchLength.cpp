#include "GibbsBranchLength.hpp"


GibbsBranchLength::GibbsBranchLength(std::unique_ptr<LikelihoodEvaluator> _eval)
  : BranchLengthMultiplier(0)
  , eval(std::move(_eval))
{
  name = "estGibbsBL"; 
  category = Category::BRANCH_LENGTHS; 
  relativeWeight = 20 ; 
} 



GibbsBranchLength::GibbsBranchLength(const GibbsBranchLength& rhs) 
  : BranchLengthMultiplier(rhs)
  , eval(rhs.eval->clone())
{
}


void GibbsBranchLength::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
{
  Branch b = proposeBranch(traln, rand);     
  
  auto params = getBranchLengthsParameterView();
  assert(params.size() == 1 );     
  auto param = params[0]; 

  b  = traln.getBranch(b,param); 
  double initBl = b.getLength(param); 
  savedBranch = b; 

  GibbsProposal::drawFromEsitmatedPosterior(b, *eval, traln, rand,  MAX_ITER, hastings, param); 
  double newZ = b.getLength(param);

  traln.setBranch(b, param); 

  prior.updateBranchLengthPrior(traln,initBl, newZ, param);  
}
