#ifndef DIVRATESLIDER_H
#define DIVRATESLIDER_H


#include "AbstractProposal.hpp"

  // AbstractProposal( Category cat, std::string  _name, double weight, double minTuning, double maxTuning, bool needsFullTraversal )  ; 


class DivRateSlider : public AbstractProposal
{

public: 
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
  DivRateSlider()
    : AbstractProposal(Category::DIVERGENCE_TIMES, "divRateSlider", 10 , 1e-5, 1e2, false)
  {}

  

  virtual ~DivRateSlider(){}
};




#endif /* DIVRATESLIDER_H */
