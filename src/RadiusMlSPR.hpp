#ifndef _RADIUSMLSPR_H
#define _RADIUSMLSPR_H

#include "AbstractProposal.hpp"
#include "Path.hpp"

typedef struct  _insertWeight
{
  branch b; 
  double lnl; 
  double weightInFirst; 
  double weightInSecond; 
  boolean containedInFirst; 
  boolean containedInSecond; 
  double ratio; 		/* TODO ? */
  struct _insertWeight *next; 
} insertList ; 



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
  Path path; 

  int radius; 
  int ratio; 
  
}; 



#endif
