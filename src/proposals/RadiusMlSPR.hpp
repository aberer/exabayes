#ifndef _RADIUSMLSPR_H
#define _RADIUSMLSPR_H


#include "AbstractProposal.hpp"
#include "Path.hpp"

typedef struct  _insertWeight
{
  Branch b; 
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
  RadiusMlSPR( int radius);
  virtual ~RadiusMlSPR(){}
  
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(  LikelihoodEvaluatorPtr &evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  

private: 
  Path path; 

  int radius; 
  int ratio; 
  
}; 



#endif

