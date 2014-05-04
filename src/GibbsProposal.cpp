#include "GibbsProposal.hpp"
#include "proposals/AbstractProposal.hpp"
#include "priors/ExponentialPrior.hpp"
#include "system/BoundsChecker.hpp"
#include "math/Density.hpp"
#include "AdHocIntegrator.hpp"
#include "math/Arithmetics.hpp"

// between .85 and 1 
// could be between 300 and 500 

#define VERBOSE_INFO

const double GibbsProposal::EST_FAC = 1.5; 


BranchLength GibbsProposal::drawFromEsitmatedPosterior(const BranchLength &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, log_double &hastings,  const AbstractParameter*  blParam) 
{
  auto result = optimiseBranch(traln, branch, eval, maxIter ,blParam);
  auto optLen = result[0]; 
  auto nrD1  = result[1]; 
  auto nrD2 = result[2]; 
  auto resultBranch = branch; 
  resultBranch.setLength(optLen); 

  auto proposalResult = propose(resultBranch.getInterpretedLength(traln, blParam), nrD1, nrD2, rand); 
  auto proposal = proposalResult[0];  
  auto alpha = proposalResult[1]; 
  auto beta = proposalResult[2]; 

  resultBranch.setConvertedInternalLength(traln,blParam, proposal);

  if(not BoundsChecker::checkBranch(resultBranch))
    BoundsChecker::correctBranch(resultBranch); 

  auto lnBackP = Density::lnGamma(branch.getInterpretedLength(traln, blParam), alpha, beta); 
  auto lnForP = Density::lnGamma(resultBranch.getInterpretedLength(traln, blParam), alpha, beta);

  // assert(lnBackP != 0 && lnForP != 0 ); 

  hastings *= lnBackP / lnForP; 

  return resultBranch; 
}



std::pair<double,double> GibbsProposal::getGammaParams(double nrOpt, double nrD1, double nrD2)
{
#ifdef VERBOSE_INFO
  // tout << MAX_SCI_PRECISION << nrOpt << "\t" << nrD2 << std::endl; 
#endif

  bool hasConverged = not (nrOpt < 1e-5 && 1e-1 < fabs(nrD1)) ; 

  double alpha = 0,  beta = 0; 
  if(not hasConverged)
    {
      alpha = 1 ; 
      beta = nrD1 ; 

    }
  else 
    {
      double c = EST_FAC  /  (- nrD2) ; 
      beta =  ((nrOpt + sqrt(nrOpt * nrOpt + 4 * c)) / (2 * c )) ; 
      alpha = nrOpt * beta + 1 ; 
    }

  assert(nrD2 < 0 || not hasConverged ); 

  return make_pair(alpha, beta);
}


std::array<double,3> GibbsProposal::propose(double nrOpt, double nrd1, double nrd2, Randomness &rand)
{
  auto params = getGammaParams(nrOpt, nrd1, nrd2); 
  double proposal = rand.drawRandGamma(params.first, params.second); 
  
  auto result = std::array<double,3>{{ proposal, params.first, params.second} }; 
  return result; 
}


std::array<double,3> GibbsProposal::optimiseBranch( TreeAln &traln, const BranchLength& b, LikelihoodEvaluator& eval, int maxIter,  AbstractParameter const  *param)
{
  double nrD1 = 0; 
  double nrD2 = 0; 

#if 0
  auto p = b.findNodePtr(traln ) ,
    q = p->back ; 
#endif

  assert(traln.getNumberOfPartitions() == 1 ); 

#if 0 
  eval.evalSubtree(traln, 0, b.toPlain() ); 
  eval.evalSubtree(traln, 0, b.toPlain().getInverted());
#endif
  
  auto prior = param->getPrior(); 
  assert(dynamic_cast<ExponentialPrior*> (prior) != nullptr); 
  double result = 0;   

#if 0 

  double lambda = dynamic_cast<ExponentialPrior*>(param->getPrior())->getLamda(); 

  double init = b.getLength(); 


  double add2nrd1 = lambda * traln.getMeanSubstitutionRate( { 0 }); 

  assert(0); 

  makenewzGeneric(&(traln.getTrHandle()), &(traln.getPartitionsHandle()), p, q, &init, maxIter,  &result , &nrD1,  &nrD2, add2nrd1, PLL_FALSE, &eval.getArrayReservoir()) ;

#endif

  return std::array<double,3>{{result, nrD1, nrD2}}; 
}
