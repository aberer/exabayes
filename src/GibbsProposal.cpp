#include "GibbsProposal.hpp"
#include "AbstractProposal.hpp"
#include "densities.h"

const double GibbsProposal::EST_EXP = -1; 

// const double GibbsProposal::EST_FAC = 0.015; 
const double GibbsProposal::EST_FAC = 1; 



/** 
    @brief draws a branch from the estimated posterior probability   
    
    internal length used as a starting point for newton-raphson
    optimization

    important: uses the current length for the hastings

    @param Branch branch -- contains the initial branch length
    
 */ 
void GibbsProposal::drawFromEsitmatedPosterior(Branch &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, double &hastings) 
{
  double initLength = branch.getLength(); 

  // std::cout <<  << std::endl; 
  
  double nrD2 = 0; 
  optimiseBranch(traln, branch, eval, nrD2, maxIter);

  if( nrD2 > 0 )
    {
      std::cerr << "2nd derivate > 0. HACK: setting it to negative value. Resolve this!" << std::endl; 
      nrD2 = - nrD2; 
    }

  double nrOpt = branch.getInterpretedLength(traln); 
  double c = EST_FAC  / -nrD2 ; 
  double beta =  ((nrOpt + sqrt(nrOpt * nrOpt + 4 * c)) / (2 * c )); 
  double alpha = nrOpt * beta + 1 ; 

  double proposal = rand.drawRandGamma(alpha,   beta);
  double prop = branch.getInternalLength(traln, proposal );   

  branch.setLength(initLength); 
  double lnBackP = logGammaDensity(branch.getInterpretedLength(traln), alpha, beta); 
  double lnForP = logGammaDensity(proposal, alpha, beta); 

  // std::cout << std::endl <<  "initL=" << branch.getInterpretedLength(traln) << "\tnrOpt=" << nrOpt << "\tnrD2=" << nrD2 << "\tproposal=" 
  // 	    << proposal <<  "\thastings=" << lnBackP - lnForP   <<std::endl; 

  
  // hastings += (lnBackP - lnForP); 
  AbstractProposal::updateHastings(hastings,   exp(lnBackP - lnForP) , "theGibbs");
  branch.setLength(prop); 
}



void GibbsProposal::optimiseBranch( TreeAln &traln, Branch &b, LikelihoodEvaluator& eval, double &secDerivative, int maxIter)
{
  auto p = b.findNodePtr(traln ), 
    q = p->back; 

  if(not p->x == 1 )
    eval.evalSubtree(traln, b); 
  if(not q->x == 1 )
    eval.evalSubtree(traln,b.getInverted());

  // eval.evaluate(traln,b, false);


  double lambda = 10;
  double result = 0;   
  double init = b.getLength(); 
  double firstDerivative = 0 ; 

#if HAVE_PLL != 0
  makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#else 
  makenewzGeneric(traln.getTr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#endif

  b.setLength(result); 
  
  // std::cout<< "result=" << result << "\tresult_iter=" << b.getInterpretedLength(traln) << "\tfirst=" << firstDerivative<< "\tsecond="<< secDerivative << "\tafter " << maxIter << std::endl; 
}
