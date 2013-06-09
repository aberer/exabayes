/** 
    @file PartitionProposal.hpp
    
    @brief proposes to change ONE parameter (e.g., revmat, shape, stateFreq)
    
    This parameter may be linked accross partitions.     
 */ 


#ifndef _PART_PROPO_H
#define _PART_PROPO_H

#include "AbstractProposal.hpp"
#include "eval.h"
#include "tune.h"

template<typename FUN, typename PARAM>
class PartitionProposal : public AbstractProposal
{
public: 
  PartitionProposal(  double _param, string _name); 
  virtual ~PartitionProposal(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 

  virtual void autotune();

  virtual PartitionProposal<FUN,PARAM>* clone() const; 

  static double relativeWeight;

  virtual double getRelativeWeight() const {return relativeWeight; }

private: 
  double parameter; 		
  vector<double> values; 

  vector<branch> storedBranches; 

}; 

#include "PartitionProposalImpl.hpp"

#endif
