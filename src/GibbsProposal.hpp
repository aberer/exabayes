#ifndef _GIBBS_PROPOSAL 
#define _GIBBS_PROPOSAL 

#include <array>
#include "LikelihoodEvaluator.hpp"


class GibbsProposal
{
public: 
  static std::array<double,3> propose(double nropt, double nrd2, Randomness &rand); 
  /** 
      @brief draw a new branch length from the posterior branch length distribution of a branch 
   */ 
  static BranchLength drawFromEsitmatedPosterior(const BranchLength &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, double &hastings,  const AbstractParameter* param) ; 
  /** 
      @brief optimises the branch 
      @return the branch with optimised branch lengths for all parameters 
   */ 
  // static BranchLength optimiseBranch( TreeAln &traln, const Branch& b, LikelihoodEvaluator& eval,  int maxIter,  const std::vector<AbstractParameter*> &params); 
  /** 
      @brief optimises the branch
      @return the optimised branch 
   */ 
  static BranchLength optimiseBranch( TreeAln &traln, const BranchLength& b, LikelihoodEvaluator& eval, double &firstDerivative, double &secDerivative, int maxIter,  const AbstractParameter*  param); 

private: 
  const static double EST_EXP; 
  const static double EST_FAC; 

}; 

#endif
