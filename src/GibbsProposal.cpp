#include "GibbsProposal.hpp"
#include "AbstractProposal.hpp"
#include "priors/ExponentialPrior.hpp"
#include "BoundsChecker.hpp"
#include "Density.hpp"
#include "AdHocIntegrator.hpp"
#include "Arithmetics.hpp"

// between .85 and 1 
#define FIXED_SMALL_ALPHA  0.95
// could be between 300 and 500 
#define FIXED_SMALL_BETA  450

// const double GibbsProposal::EST_FAC = 0.015; 
// const double GibbsProposal::EST_FAC = 1.5; 

// const double GibbsProposal::EST_EXP = -1; 
const double GibbsProposal::EST_FAC = 1.5; 


BranchLength GibbsProposal::drawFromEsitmatedPosterior(const BranchLength &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, double &hastings,  const AbstractParameter*  blParam) 
{
  auto result = optimiseBranch(traln, branch, eval, maxIter ,blParam);
  auto optLen = result[0]; 
  auto nrD2 = result[2]; 
  auto resultBranch = branch; 
  resultBranch.setLength(optLen); 

  auto proposalResult = propose(resultBranch.getInterpretedLength(traln, blParam), nrD2, rand ); 
  auto proposal = proposalResult[0];  
  auto alpha = proposalResult[1]; 
  auto beta = proposalResult[2]; 

  resultBranch.setConvertedInternalLength(traln,blParam, proposal);

  if(not BoundsChecker::checkBranch(resultBranch))
    BoundsChecker::correctBranch(resultBranch); 


  double lnBackP = Density::lnGamma(branch.getInterpretedLength(traln, blParam), alpha, beta); 
  double lnForP = Density::lnGamma(proposal, alpha, beta); 

  assert(lnBackP != 0 && lnForP != 0 ); 

  AbstractProposal::updateHastingsLog(hastings,   lnBackP - lnForP , "theGibbs");
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
      // tout << "GIBBS-BAD" << std::endl;
    }
  else 
    {
      double c = EST_FAC  / nrd2 ; 
      beta =  ((nrOpt + sqrt(nrOpt * nrOpt + 4 * c)) / (2 * c )) ; 
      alpha = nrOpt * beta + 1 ; 

      proposal = rand.drawRandGamma(alpha,   beta);
    }

  auto result = std::array<double,3>{{ proposal, alpha, beta} }; 

  return result; 
}


std::array<double,3> GibbsProposal::optimiseBranch( TreeAln &traln, const BranchLength& b, LikelihoodEvaluator& eval, int maxIter,  AbstractParameter const  *param)
{
  double nrD1 = 0; 
  double nrD2 = 0; 

  auto p = b.findNodePtr(traln ), 
    q = p->back; 

  assert(traln.getNumberOfPartitions() == 1 ); 

  eval.evalSubtree(traln, 0, b.toPlain() ); 
  eval.evalSubtree(traln, 0, b.toPlain().getInverted());
  
  auto prior = param->getPrior(); 
  assert(dynamic_cast<ExponentialPrior*> (prior) != nullptr); 
  double lambda = dynamic_cast<ExponentialPrior*>(param->getPrior())->getLamda(); 

  double result = 0;   
  double init = b.getLength(); 

#if HAVE_PLL != 0
  makenewzGeneric(traln.getTr(), traln.getPartitionsPtr(), p, q, &init, maxIter,  &result , &nrD1,  &nrD2, lambda, FALSE) ;
#else 
  makenewzGeneric(traln.getTr(), p, q, &init, maxIter,  &result , &nrD1,  &nrD2, lambda, FALSE) ;
#endif

  return std::array<double,3>{{result, nrD1, nrD2}}; 
}
