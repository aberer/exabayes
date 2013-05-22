#ifndef _EXTENDED_SPR_H
#define _EXTENDED_SPR_H

#include "axml.h"
#include "AbstractProposal.hpp"
#include "Randomness.hpp"
#include "Path.hpp"
#include "SprMove.hpp"


class ExtendedSPR : public AbstractProposal
{
public: 
  ExtendedSPR(  double stopProb, double multiplier); 
  virtual ~ExtendedSPR(); 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune() {}	// disabled 
  virtual AbstractProposal* clone() const;  

  static double relativeWeight;

  virtual double getRelativeWeight() const {return relativeWeight; }

protected: 
  double stopProb; 
  double multiplier; 
  Path modifiedPath; 

  SprMove move; 

  void drawPathForESPR( TreeAln& traln, Randomness &rand, double stopProp ); 
}; 


#endif
