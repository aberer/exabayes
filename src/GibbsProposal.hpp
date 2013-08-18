#ifndef _GIBBS_PROPOSAL 
#define _GIBBS_PROPOSAL 

#include "LikelihoodEvaluator.hpp"


class GibbsProposal
{
public: 
  static Branch drawFromEsitmatedPosterior(const Branch &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, double &hastings,  AbstractParameter* const  param) ; 

private: 
  static Branch optimiseBranch( TreeAln &traln, Branch b, LikelihoodEvaluator& eval, double &firstDerivative, double &secDerivative, int maxIter,  AbstractParameter* const  param); 

  const static double EST_EXP; 
  const static double EST_FAC; 

}; 

#endif
