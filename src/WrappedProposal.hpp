
#include "proposalFunction.h"

#include "AbstractProposal.hpp"

class Chain; 



class WrappedProposal : public AbstractProposal
{
public: 
  WrappedProposal(proposalFunction* pf, Chain *chain) ;
  virtual ~WrappedProposal(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune() ;
  virtual WrappedProposal* clone() const; 

  virtual void setOwningChain(Chain *_chain)  {chain = _chain; }

private: 
  proposalFunction *pfun; 
  Chain *chain; 

};
