
#include "bayes.h"

#include "AbstractProposal.hpp"

class Chain; 



class WrappedProposal : public AbstractProposal
{
public: 
  WrappedProposal(proposalFunction* pf, Chain *chain) ;
  virtual ~WrappedProposal(){}

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) ; 
  virtual void resetState(TreeAln &traln, PriorManager &prior) ; 
  virtual void autotune() ;

  virtual void setOwningChain(Chain *_chain)  {chain = _chain; }

private: 
  proposalFunction *pfun; 
  Chain *chain; 

};
