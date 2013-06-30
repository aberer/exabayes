#include "GibbsProposal.hpp"
#include "AbstractProposal.hpp"
#include "densities.h"

const double GibbsProposal::EST_EXP = -1; 
const double GibbsProposal::EST_FAC = 0.015; 


/** 
    @brief draws a branch from the estimated posterior probability   
    
    internal length used as a starting point for newton-raphson
    optimization

    important: uses the current length for the hastings

    @param Branch branch -- contains the initial branch length
    
 */ 
void GibbsProposal::drawFromEsitmatedPosterior(Branch &branch, LikelihoodEvaluatorPtr& eval, TreeAln &traln, Randomness& rand, double initVal,  int maxIter, double &hastings) 
{
  auto p = branch.findNodePtr(traln); 
  double initLength = branch.getLength(); 
  
  double nrD2 = 0; 
  optimiseBranch(traln, branch, eval, nrD2, maxIter);

  double nrOpt = branch.getInterpretedLength(traln); 
  double c = EST_FAC  / -nrD2 ; 
  double beta =  ((nrOpt + sqrt(nrOpt * nrOpt + 4 * c)) / (2 * c )); 
  double alpha = nrOpt * beta + 1 ; 

  double proposal = rand.drawRandGamma(alpha,   beta);
  double prop = branch.getInternalLength(traln, proposal );   

  branch.setLength(initLength); 
  double lnBackP = logGammaProb(branch.getInterpretedLength(traln), alpha, beta); 
  double lnForP = logGammaProb(proposal, alpha, beta); 

  AbstractProposal::updateHastings(hastings,   exp(lnBackP - lnForP) , "theGibbs");
  branch.setLength(prop); 
}



void GibbsProposal::optimiseBranch( TreeAln &traln, Branch &b, LikelihoodEvaluatorPtr &eval, double &secDerivative, int maxIter)  
{
  auto p = b.findNodePtr(traln ), 
    q = p->back; 

  if(not p->x == 1 )
    eval->evalSubtree(traln, b); 
  if(not q->x == 1 )
    eval->evalSubtree(traln,b.getInverted());
    
    eval->evaluate(traln,b, false);

#ifdef TODO
  // assert(0); 
#endif

  double lambda = 10;
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
