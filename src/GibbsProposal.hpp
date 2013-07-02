#ifndef _GIBBS_PROPOSAL 
#define _GIBBS_PROPOSAL 

#include "LikelihoodEvaluator.hpp"


class GibbsProposal
{
public: 
  static void drawFromEsitmatedPosterior(Branch &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand, double initVal,  int maxIter, double &hastings) ; 

private: 
  static void optimiseBranch( TreeAln &traln, Branch &b, LikelihoodEvaluator &eval, double &secDerivative, int maxIter)  ; 

  const static double EST_EXP; 
  const static double EST_FAC; 

}; 

#endif
