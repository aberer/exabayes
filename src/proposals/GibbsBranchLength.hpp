#include "GibbsProposal.hpp" 
#include "BranchLengthMultiplier.hpp"
#include "Category.hpp"


#define MAX_ITER 30 

class GibbsBranchLength : public BranchLengthMultiplier
{
public: 
  GibbsBranchLength(std::unique_ptr<LikelihoodEvaluator> _eval); 
  GibbsBranchLength(const GibbsBranchLength& rhs); 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void autotune(){}

  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { assert(0); return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}

  virtual AbstractProposal* clone() const  { return new GibbsBranchLength(*this); }  

private: 
  std::unique_ptr<LikelihoodEvaluator> eval; 
}; 
