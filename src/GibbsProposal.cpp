#include "GibbsProposal.hpp"
#include "AbstractProposal.hpp"
#include "priors/ExponentialPrior.hpp"
#include "BoundsChecker.hpp"
#include "densities.h"
#include "AdHocIntegrator.hpp"
#include "Arithmetics.hpp"

// between .85 and 1 
#define FIXED_SMALL_ALPHA  0.95
// could be between 300 and 500 
#define FIXED_SMALL_BETA  450

// const double GibbsProposal::EST_FAC = 0.015; 
// const double GibbsProposal::EST_FAC = 1.5; 

const double GibbsProposal::EST_EXP = -1; 
const double GibbsProposal::EST_FAC = 1.5; 


/** 
    @brief draws a branch from the estimated posterior probability   
    
    internal length used as a starting point for newton-raphson
    optimization

    important: uses the current length for the hastings

    @param Branch branch -- contains the initial branch length
    
 */ 
BranchLength GibbsProposal::drawFromEsitmatedPosterior(const BranchLength &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, double &hastings,  const AbstractParameter*  blParam) 
{
  double nrD2 = 0; 
  double nrD1 = 0; 
  auto optimizedBranch = optimiseBranch(traln, branch, eval, nrD1, nrD2, maxIter ,blParam);

  auto proposalResult = propose(optimizedBranch.getInterpretedLength(traln, blParam), nrD2, rand ); 
  auto proposal = proposalResult[0];  
  auto alpha = proposalResult[1]; 
  auto beta = proposalResult[2]; 

  auto resultBranch = optimizedBranch; 
  resultBranch.setConvertedInternalLength(traln,blParam, proposal);
  
// #ifdef _EXPERIMENTAL_INTEGRATION_MODE
  // ahInt->prepareForBranch(branch, traln); 
  // auto samples = ahInt->integrate(branch, traln, 1000, 10); 
  // auto meanVar = Arithmetics::getMeanAndVar(samples); 
  // double var = meanVar.second; 

  // tout << "INT\t"  << MAX_SCI_PRECISION << nrOpt << "\t" 
  //      << nrD2 << "\t"
  //      <<  var; 
// #endif

  if(not BoundsChecker::checkBranch(resultBranch))
    BoundsChecker::correctBranch(resultBranch); 


  double lnBackP = logGammaDensity(branch.getInterpretedLength(traln, blParam), alpha, beta); 
  double lnForP = logGammaDensity(proposal, alpha, beta); 

  assert(lnBackP != 0 && lnForP != 0 ); 

  AbstractProposal::updateHastings(hastings,   exp(lnBackP - lnForP) , "theGibbs");
  return resultBranch; 
}


std::array<double,3> GibbsProposal::propose(double nrOpt, double nrd2, Randomness &rand)
{
  bool d2isPositive = nrd2 > 0 ; 
  if(not d2isPositive) 
    nrd2 = -nrd2; 

  double proposal = 0; 
  double alpha = 0,  beta = 0; 
  if(d2isPositive)
    {
      alpha = FIXED_SMALL_ALPHA; 
      beta = FIXED_SMALL_BETA; 
      proposal = rand.drawRandGamma(alpha, beta); 
      tout << "GIBBS-BAD" << std::endl;
    }
  else 
    {
      // nrOpt = optimizedBranch.getInterpretedLength(traln, blParam); 
      double c = EST_FAC  / nrd2 ; 
      beta =  ((nrOpt + sqrt(nrOpt * nrOpt + 4 * c)) / (2 * c )) ; 
      alpha = nrOpt * beta + 1 ; 

      proposal = rand.drawRandGamma(alpha,   beta);
    }

  return std::array<double,3>{{proposal, alpha, beta}}; 
}


BranchLength GibbsProposal::optimiseBranch( TreeAln &traln, const BranchLength& b, LikelihoodEvaluator& eval, double &firstDerivative, double &secDerivative, int maxIter,  AbstractParameter const  *param)
{
  auto resultBranch = b; 
  auto p = b.findNodePtr(traln ), 
    q = p->back; 

  if( p->x != 1 )
    eval.evalSubtree(traln, b.toPlain()); 
  if( q->x != 1 )
    eval.evalSubtree(traln,b.toPlain().getInverted());

  // eval.evaluate(traln,b, false);

  auto prior = param->getPrior(); 
  assert(dynamic_cast<ExponentialPrior*> (prior) != nullptr); 
  double lambda = dynamic_cast<ExponentialPrior*>(param->getPrior())->getLamda(); 

  double result = 0;   
  double init = b.getLength(); 
  // double firstDerivative = 0 ; 

#if HAVE_PLL != 0
  makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#else 
  makenewzGeneric(traln.getTr(), p, q, &init, maxIter,  &result , &firstDerivative,  &secDerivative, lambda, FALSE) ;
#endif

  resultBranch.setLength(result); 
  
  return resultBranch; 
  // std::cout<< "result=" << result << "\tresult_iter=" << b.getInterpretedLength(traln) << "\tfirst=" << firstDerivative<< "\tsecond="<< secDerivative << "\tafter " << maxIter << std::endl; 
}


// BranchLength GibbsProposal::optimiseBranch( TreeAln &traln, const BranchLength& b, LikelihoodEvaluator& eval,  int maxIter,  const std::vector<AbstractParameter*> &params)
// {
//   auto result = b; 
  
//   for(auto param : params)
//     {
//       double a,c; 
//       auto partResult = optimiseBranch(traln, b, eval, a,c, maxIter, param); 
//       result.setLength(partResult.getLength(param),param);
//     }

//   return result; 
// }




