/** 
    @file ParameterProposal.hpp
    
    @brief represents proposals that change the alpha parameter of the
    gamma distribution of rate heterogeneity.
    
 */ 


#ifndef _SIMPLERATEHET_H
#define _SIMPLERATEHET_H

#include "AbstractProposal.hpp"

template<typename FUN, typename PARAM>
class PartitionProposal : public AbstractProposal
{
public: 
  PartitionProposal( double relativeWeight, double _param, string _name); 
  virtual ~PartitionProposal(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune();

  virtual PartitionProposal<FUN,PARAM>* clone() const; 

private: 
  int model; 			// which model
  double parameter; 
  vector<double> values; 
}; 

#include "PartitionProposalImpl.hpp"

#endif
