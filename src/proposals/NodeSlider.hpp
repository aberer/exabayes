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
  virtual void evaluateProposal(  LikelihoodEvaluator *evaluator, TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  
  Branch proposeBranch(const TreeAln &traln, Randomness &rand) const ; 
  Branch proposeOtherBranch(const Branch &firstBranch, const TreeAln& traln, Randomness& rand) const ; 
  
  // or have a generic set-evaluation method? 
  void prepareForEvaluation( TreeAln &traln) const ; 
  std::pair<Branch,Branch> prepareForSetExecution(TreeAln& traln, Randomness &rand) ; 
  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  

  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out) const { } //disabled

protected: 
  double multiplier; 
  Branch oneBranch; 
  Branch otherBranch; 
};  

#endif
