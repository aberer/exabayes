#include "bayes.h"
#include "axml.h"

#include "AbstractProposal.hpp"

class Chain; 


class ExtendedTBR : public AbstractProposal
{
public: 
  ExtendedTBR( double relativeProbability, double extensionProb);
  virtual ~ExtendedTBR(){};

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand); 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior); 
  virtual void resetState(TreeAln& traln, PriorManager &prior); 
  virtual void autotune();
  
  virtual void setOwningChain(Chain *chain){}
  
private: 
  double extensionProbability; 


}; 


