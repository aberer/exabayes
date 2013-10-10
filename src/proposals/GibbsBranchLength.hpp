#include "BranchLengthMultiplier.hpp"


// #define _EXPERIMENTAL

#define MAX_ITER 30 

class GibbsBranchLength : public BranchLengthMultiplier
{
public: 
  GibbsBranchLength(bool doTwo);

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
#ifdef  _EXPERIMENTAL
  virtual void resetState(TreeAln &traln) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, const BranchPlain &branchSuggestion) ; 
#endif
  
  virtual void autotune(){}

  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { assert(0); return std::make_pair(BranchPlain(0,0),BranchPlain(0,0) );}

  virtual AbstractProposal* clone() const  { return new GibbsBranchLength(*this); }  

private: 
  BranchLength extraBranch; 
  bool doTwo ;  
}; 
