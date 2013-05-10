#ifndef _NODESLIDER_H
#define _NODESLIDER_H

#include "AbstractProposal.hpp" 
#include "Path.hpp"


class NodeSlider : public AbstractProposal
{
public:   
  NodeSlider(Chain *chain, double relativeProbability, double multiplier);
  virtual ~NodeSlider(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune() {}	// disabled 
  virtual void setOwningChain(Chain *_chain) {chain = _chain;}

  virtual AbstractProposal* clone() const;  

protected: 
  Chain* chain; 
  double multiplier; 

  Path path; 

};  

#endif
