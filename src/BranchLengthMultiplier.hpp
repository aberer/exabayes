#ifndef _BRANCHLENGTHSMULTIPLIER_H
#define _BRANCHLENGTHSMULTIPLIER_H

#include "AbstractProposal.hpp"



class BranchLengthMultiplier : public AbstractProposal
{
public: 
  BranchLengthMultiplier(Chain *chain, double relativeWeight, double multiplier); 
  virtual ~BranchLengthMultiplier(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune();

  virtual void setOwningChain(Chain *_chain) {chain = _chain;}

  virtual AbstractProposal* clone() const;  
  
private: 
  Chain *chain; 
  double multiplier;  
  
  double savedZ;  
}; 


#endif
