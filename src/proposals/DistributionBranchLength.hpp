#ifndef _DISTRBUTION_BRANCH_LENGTH
#define  _DISTRBUTION_BRANCH_LENGTH

#include "BranchLengthMultiplier.hpp"

#include "DistributionProposer.hpp"
#include "GammaProposer.hpp"

#include "ComplexTuner.hpp"

#include "OptimizedParameter.hpp"

template <class C>
class DistributionBranchLength : public BranchLengthMultiplier
{
public: 
  DistributionBranchLength();

  virtual void createProposer(TreeAln &traln, LikelihoodEvaluator& eval, BranchPlain b); 
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void autotune();  
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand) ; 
  virtual AbstractProposal* clone() const  { return new DistributionBranchLength(*this); }  
  virtual void extractProposer(TreeAln& traln, const OptimizedParameter& param) 
  {
    _proposer = param.getProposerDistribution< C >(traln, _convTuner.getParameter(), _nonConvTuner.getParameter() );
  }

  virtual void accept(); 
  virtual void reject(); 

private: 
  DistributionProposer<C> _proposer;

  bool _converged; 

  ComplexTuner _convTuner; 
  ComplexTuner _nonConvTuner; 
}; 



void calcDiagptable(const double z, const int states, const int numberOfCategories, const double *rptr, const double *EIGN, double *diagptable); 
#include "proposals/DistributionBranchLength.tpp"

#endif
