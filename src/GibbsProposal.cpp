#include "GibbsProposal.hpp"

#include "AbstractProposal.hpp"

const double GibbsProposal::EST_EXP = -0.87246; 
const double GibbsProposal::EST_FAC = 0.008067; 


/** 
    @brief draws a branch from the estimated posterior probability   
    
    internal length used as a starting point for newton-raphson
    optimization

    important: uses the current length for the hastings

    @param Branch branch -- contains the initial branch length
    
 */ 
double GibbsProposal::drawFromEsitmatedPosterior(Branch &branch, LikelihoodEvaluatorPtr& eval, TreeAln &traln, Randomness& rand, double initVal,  int maxIter, double &hastings) 
{
  auto p = branch.findNodePtr(traln); 
  double initLength = branch.getLength(); 

  
  
  double nrD2 = 0; 
  optimiseBranch(traln, branch, eval, nrD2, maxIter);

  double nrOpt = branch.getLength(); 
  double c = EST_FAC * pow(-nrD2, EST_EXP); 
  double beta = nrOpt + sqrt(pow(nrOpt,2) + 4 * c) / (2 * c ); 
  double alpha = nrOpt  * beta + 1 ; 
  
  double proposal = rand.drawRandGamma(alpha, beta);
  double prop = branch.getInternalLength(traln, proposal );   
  AbstractProposal::updateHastings(hastings, log(initLength) / log(prop), "theGibbs");
  branch.setLength(prop); 
  
  return 0 ; 
}



void GibbsProposal::optimiseBranch( TreeAln &traln, Branch &b, LikelihoodEvaluatorPtr &eval, double &secDerivative, int maxIter)  
{
  auto p = b.findNodePtr(traln ), 
    q = p->back; 

  if(not p->x == 1 )
    eval->evalSubtree(traln, b); 
  if(not q->x == 1 )
    eval->evalSubtree(traln,b.getInverted());
  
#ifdef TODO
  // assert(0); 
#endif

  double lambda = 10; 		// unnecesary ? TODO   
  double result = 0;   
  double init = b.getLength(); 
  double firstDerivative = 0 ; 

#if HAVE_PLL != 0
    makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#else 
  assert(0); 
#endif

  b.setLength(result); 
}
