#ifndef _NODESLIDER_H
#define _NODESLIDER_H

#include "AbstractProposal.hpp" 
#include "Path.hpp"


class NodeSlider : public AbstractProposal
{
public:   
  NodeSlider( double multiplier);
  virtual ~NodeSlider(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  

  virtual void readFromCheckpointCore(std::ifstream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ofstream &out) { } //disabled

protected: 
  double multiplier; 
  // Path path; 
  Branch oneBranch; 
  Branch otherBranch; 

};  

#endif
