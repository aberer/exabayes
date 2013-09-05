#ifndef _NODESLIDER_H
#define _NODESLIDER_H

#include "AbstractProposal.hpp" 
#include "Path.hpp"


class NodeSlider : public AbstractProposal
{
public:   
  NodeSlider( double multiplier);
  virtual ~NodeSlider(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) ; 
  virtual void resetState(TreeAln &traln) ; 
  
  BranchPlain proposeBranch(const TreeAln &traln, Randomness &rand) const ; 
  BranchPlain proposeOtherBranch(const BranchPlain &firstBranch, const TreeAln& traln, Randomness& rand) const ; 
  
  // or have a generic set-evaluation method? 
  void prepareForEvaluation( TreeAln &traln) const ; 
  std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln& traln, Randomness &rand) ; 
  virtual void autotune() {}	// disabled 

  virtual AbstractProposal* clone() const;  

  virtual void readFromCheckpointCore(std::istream &in) {   } // disabled
  virtual void writeToCheckpointCore(std::ostream &out) const { } //disabled

protected: 
  double multiplier; 
  BranchLength oneBranch; 
  BranchLength otherBranch; 
};  

#endif
