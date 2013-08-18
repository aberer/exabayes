#ifndef _BRANCHLENGTHSMULTIPLIER_H
#define _BRANCHLENGTHSMULTIPLIER_H

#include "AbstractProposal.hpp"

#include <limits>


class BranchLengthMultiplier : public AbstractProposal
{
public: 
  BranchLengthMultiplier(  double multiplier); 
  virtual ~BranchLengthMultiplier(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln) ; 
  virtual void resetState(TreeAln &traln) ; 

  virtual void autotune();

  virtual AbstractProposal* clone() const;  

  // virtual Branch prepareForSetExecution(TreeAln &traln, Randomness &rand)  ;
  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand) ;

  virtual Branch proposeBranch(const TreeAln &traln, Randomness &rand) const ;   

  virtual void readFromCheckpointCore(std::istream &in); 
  virtual void writeToCheckpointCore(std::ostream &out) const; 

protected: 
  double multiplier;  
  Branch savedBranch;   

}; 


#endif
