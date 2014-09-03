#ifndef DIVTIMESLIDER_H
#define DIVTIMESLIDER_H


#include "AbstractProposal.hpp"

class DivTimeSlider : public AbstractProposal
{ 
public:				// INHERITED 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) ; 
  virtual void resetState(TreeAln &traln)  ; 
  virtual void autotune()   ;
  virtual AbstractProposal* clone() const ;  
  virtual BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const ; 
  virtual std::vector<nat> getInvalidatedNodes(const TreeAln &traln) const ;  
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  ; 
  virtual void writeToCheckpointCore(std::ostream &out)const   ;  
  virtual void readFromCheckpointCore(std::istream &in) ; 

public:
  DivTimeSlider()
    : AbstractProposal(Category::DIVERGENCE_TIMES, "divTimeSlider", 10 , 1e-5, 1e2, false)
  {
  }

  virtual ~DivTimeSlider(){}

};




#endif /* DIVTIMESLIDER_H */
