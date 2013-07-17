#ifndef _MYTEMPLATEPROPOSAL_H  	// use whatever you want as an include guard (must NOT occur twice! )
#define _MYTEMPLATEPROPOSAL_H

#include "AbstractProposal.hpp"


class BranchCollapser : public AbstractProposal
{
public: 
  BranchCollapser();

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior); 
  virtual void autotune()  { }
  virtual AbstractProposal* clone() const ;  

  virtual void readFromCheckpointCore(std::ifstream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ofstream &out)  { } //disabled

private: 
  Branch modifiedBranch;   
}; 


#endif
