#ifndef _GIBBS_PROPOSAL 
#define _GIBBS_PROPOSAL 

#include <array>
#include "eval/LikelihoodEvaluator.hpp"


class GibbsProposal
{
public: 
  static std::array<double,3> propose(double nropt, double nrd1, double nrd2, Randomness &rand); 
  /** 
      @brief draw a new branch length from the posterior branch length distribution of a branch 
   */ 
  static BranchLength drawFromEsitmatedPosterior(const BranchLength &branch, LikelihoodEvaluator& eval, TreeAln &traln, Randomness& rand,  int maxIter, log_double &hastings,  const AbstractParameter* param) ; 
  static std::pair<double,double> getGammaParams(double nrOpt, double nrD1, double nrdD2); 
  /** 
      @brief optimises the branch
   */ 
  static std::array<double,3> optimiseBranch( TreeAln &traln, const BranchLength& b, LikelihoodEvaluator& eval, int maxIter,  AbstractParameter const  *param); 

private: 
  const static double EST_EXP; 
  const static double EST_FAC; 

}; 

#endif
