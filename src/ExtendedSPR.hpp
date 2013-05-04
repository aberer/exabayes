#ifndef _EXTENDED_SPR_H
#define _EXTENDED_SPR_H

#include "bayes.h"
#include "axml.h"
#include "AbstractProposal.hpp"
#include "Randomness.hpp"
#include "path.h"



class ExtendedSPR : public AbstractProposal
{
public: 
  ExtendedSPR( Chain *chain, double relativeWeight, double stopProb, double multiplier); 
  virtual ~ExtendedSPR(); 

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) ; 
  virtual void resetState(TreeAln &traln, PriorManager &prior) ; 

  virtual void autotune() {}	// disabled 
  virtual void setOwningChain(Chain *_chain) {chain = _chain;}

protected: 
  Chain* chain; 
  double stopProb; 
  double multiplier; 


  path* modifiedPath; 

}; 


#endif
