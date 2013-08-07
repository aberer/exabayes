#include "GibbsProposal.hpp"
#include "AbstractProposal.hpp"
#include "priors/ExponentialPrior.hpp"
#include "BoundsChecker.hpp"
#include "densities.h"

const double GibbsProposal::EST_EXP = -1; 

// const double GibbsProposal::EST_FAC = 0.015; 
const double GibbsProposal::EST_FAC = 1.5; 


/** 
    @brief draws a branch from the estimated posterior probability   
    
    internal length used as a starting point for newton-raphson
    optimization

    important: uses the current length for the hastings

    @param Branch branch -- contains the initial branch length
    
 */ 
void GibbsProposal::drawFromEsitmatedPosterior(Branch &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, double &hastings,  AbstractParameter const* blParam) 
{
  double initLength = branch.getLength(blParam); 
 
  double lambda = 1000 ; 

  double nrD2 = 0; 
  double nrD1 = 0; 
  optimiseBranch(traln, branch, eval, nrD1, nrD2, maxIter ,blParam);

  bool secondDeriWasPositive = false; 

  if( nrD2 > 0 )
    {
      // tout << "For branch: " << branch << std::endl; 
      // tout << "nrD1: " << nrD1 << "\tnrD2: " << nrD2 << std::endl; 

      // std::cerr << "2nd derivate > 0. HACK: setting it to negative value. Resolve this!" << std::endl; 
      nrD2 = - nrD2; 
      secondDeriWasPositive = true; 
    }

  double proposal = 0; 
  double alpha = 0,  beta = 0; 
  if(secondDeriWasPositive)
    {
      proposal = rand.drawRandExp(lambda);
      // tout << "drawn " << proposal << std::endl; 
    }
  else 
    {
      double nrOpt = branch.getInterpretedLength(traln, blParam); 
      double c = EST_FAC  / -nrD2 ; 
      beta =  ((nrOpt + sqrt(nrOpt * nrOpt + 4 * c)) / (2 * c )); 
      alpha = nrOpt * beta + 1 ; 
      proposal = rand.drawRandGamma(alpha,   beta);
    }

  branch.setConvertedInternalLength(traln,blParam, proposal);
  if(not BoundsChecker::checkBranch(branch))
    BoundsChecker::correctBranch(branch, blParam); 

  double prop = branch.getLength(blParam);  
  branch.setLength(initLength, blParam); 

  double lnBackP = 0,  lnForP = 0;   
  if(secondDeriWasPositive)
    {
      double old = branch.getInterpretedLength(traln, blParam); 
      tout << "using EXP" << std::endl; 
      // tout << "old= " << old << "\tproposal=" << proposal << std::endl ; 
      // tout << "p_exp(old;100)="  << exponentialDensity(old,lambda) 
      // 	   << "\tp_exp(proposal;100)=" << exponentialDensity(proposal,lambda) << std::endl; 
      // tout << "log(old)=" << log(exponentialDensity(old,lambda) )
      // 	   << "\tlog(proposal)=" << log(exponentialDensity(proposal,lambda) ) << std::endl; 
      lnBackP = log(exponentialDensity(old, lambda)); 
      lnForP = log(exponentialDensity(proposal, lambda));
      // tout << "=> hastings=" << lnBackP - lnForP << std::endl; 
    }
  else 
    {

      lnBackP = logGammaDensity(branch.getInterpretedLength(traln, blParam), alpha, beta); 
      lnForP = logGammaDensity(proposal, alpha, beta); 
    }
  
  


  assert(lnBackP != 0 && lnForP != 0 ); 

  AbstractProposal::updateHastings(hastings,   exp(lnBackP - lnForP) , "theGibbs");
  branch.setLength(prop, blParam); 
}



void GibbsProposal::optimiseBranch( TreeAln &traln, Branch &b, LikelihoodEvaluator& eval, double &firstDerivative, double &secDerivative, int maxIter, const AbstractParameter* param)
{
  auto p = b.findNodePtr(traln ), 
    q = p->back; 

  if(not p->x == 1 )
    eval.evalSubtree(traln, b); 
  if(not q->x == 1 )
    eval.evalSubtree(traln,b.getInverted());

  eval.evaluate(traln,b, false);

  auto prior = param->getPrior(); 
  assert(dynamic_cast<ExponentialPrior*> (prior) != nullptr); 
  double lambda = dynamic_cast<ExponentialPrior*>(param->getPrior())->getLamda(); 

  double result = 0;   
  double init = b.getLength(param); 
  // double firstDerivative = 0 ; 

#if HAVE_PLL != 0
  makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#else 
  makenewzGeneric(traln.getTr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#endif

  b.setLength(result, param); 
  
  // std::cout<< "result=" << result << "\tresult_iter=" << b.getInterpretedLength(traln) << "\tfirst=" << firstDerivative<< "\tsecond="<< secDerivative << "\tafter " << maxIter << std::endl; 
}
