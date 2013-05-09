/** 
    @brief our flavour of the statistical nearest neighbor interchange

    in constrast to mb, this proposal will always induce a topological
    change
    
    @notice we could multiply some additional branches
 */ 


#ifndef _STAT_NNI_H
#define _STAT_NNI_H

class Chain; 
#include "TreeAln.hpp"
#include "Path.hpp"
#include "AbstractProposal.hpp"

class StatNNI : public AbstractProposal
{
public: 
  StatNNI(Chain *_chain, double weigth, double multiplier);
  ~StatNNI(){}

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) ; 
  virtual void resetState(TreeAln &traln, PriorManager &prior) ; 

  virtual void autotune() {}	// disabled 
  virtual void setOwningChain(Chain *_chain) {chain = _chain;}

private:
  Chain* chain; 
  double multiplier; 
  Path path; 
  
};

#endif
