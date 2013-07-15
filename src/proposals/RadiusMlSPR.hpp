#ifndef _RADIUSMLSPR_H
#define _RADIUSMLSPR_H


#include "AbstractProposal.hpp"
#include "Path.hpp"

#define STRETCH_FACTOR 2 
#define WEIGHT_EPS 1e-3  	/* guided spr    */

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
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune() {}	// disabled 

  virtual void readFromCheckpointCore(std::ifstream &in) {   } 
  virtual void writeToCheckpointCore(std::ofstream &out)const { }  

  virtual AbstractProposal* clone() const;  

private: 
  Path path; 

  int radius; 
  int ratio; 
  
}; 



#endif

