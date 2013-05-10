#ifndef _RADIUSMLSPR_H
#define _RADIUSMLSPR_H

#include "AbstractProposal.hpp"


class RadiusMlSPR : public AbstractProposal
{
public: 
  RadiusMlSPR(Chain *chain, double relativeWeight, int radius);
  virtual ~RadiusMlSPR(){}
  
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune() {}	// disabled 
  virtual void setOwningChain(Chain *_chain) {chain = _chain;}


private: 
  Chain *chain; 

  int radius; 
  
}; 



#endif
