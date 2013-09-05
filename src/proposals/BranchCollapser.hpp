#ifndef _MYTEMPLATEPROPOSAL_H  	// use whatever you want as an include guard (must NOT occur twice! )
#define _MYTEMPLATEPROPOSAL_H

#include "AbstractProposal.hpp"


class BranchCollapser : public AbstractProposal
{
public: 
  BranchCollapser();

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator &eval) ; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator,TreeAln &traln) ; 
  virtual void resetState(TreeAln &traln); 
  virtual void autotune()  { }
  virtual AbstractProposal* clone() const ;  

  // virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::make_pair(BranchPlain(0,0),BranchPlain(0,0)); }

  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out)  const { } //disabled

private: 
  BranchPlain modifiedBranch;   
}; 


#endif
