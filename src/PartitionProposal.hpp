/** 
    @file ParameterProposal.hpp
    
    @brief represents proposals that change the alpha parameter of the
    gamma distribution of rate heterogeneity.
    
 */ 


#ifndef _SIMPLERATEHET_H
#define _SIMPLERATEHET_H

#include "AbstractProposal.hpp"


class Chain; 

template<typename FUN, typename PARAM>
class PartitionProposal : public AbstractProposal
{
public: 
  PartitionProposal(Chain *_chain, double relativeWeight, double _param, string _name); 
  virtual ~PartitionProposal(){}

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) ; 
  virtual void resetState(TreeAln &traln, PriorManager &prior) ; 

  virtual void autotune();	// disabled 
  virtual void setOwningChain(Chain *_chain) {chain = _chain;}

private: 
  int model; 			// which model
  double parameter; 
  vector<double> values; 
  Chain *chain ;
}; 


#include "PartitionProposalImpl.hpp"


#endif
