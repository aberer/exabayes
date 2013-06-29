#ifndef _GIBBS_PROPOSAL 
#define _GIBBS_PROPOSAL 

#include "LikelihoodEvaluator.hpp"


class GibbsProposal
{
public: 
  static double drawFromEsitmatedPosterior(Branch &branch, TreeAln &traln, Randomness& rand, double initVal,  int maxIter, double &hastings) ; 

private: 
  const static double EST_EXP; 
  const static double EST_FAC; 

}; 

#endif
