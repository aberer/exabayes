#include "AbstractProposal.hpp"

class Chain ; 


class TreeLengthMultiplier : public AbstractProposal
{
public: 
TreeLengthMultiplier(Chain *_chain, double _relativeWeight, double _multiplier); 
  ~TreeLengthMultiplier(){}

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) ; 
  virtual void resetState(TreeAln &traln, PriorManager &prior)  ; 
  virtual void autotune() ;

  virtual void setOwningChain(Chain *_chain) {chain = _chain;}


private: 
Chain *chain; 

  double multiplier; 		// the tuning variable
  
  double rememMultiplier; 	// for resetting 

void multiplyBranchLengthsRecursively(TreeAln& traln, nodeptr p, double multiHere); 


} ; 
