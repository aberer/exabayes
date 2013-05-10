/** 
    @file ParsimonySPR.hpp
    
    @brief implements a parsimony-biased SPR move similar to MrBayes. 

    @notice: MrBayes also reweights site patters -- currently we do
    not do that.
 */ 



// #include "bayes.h"
#include "axml.h"
#include "AbstractProposal.hpp"
#include "Path.hpp"

class Chain; 


class ParsimonySPR : public AbstractProposal
{
public: 
  ParsimonySPR(Chain *chain, double relativeWeight, double parsWarp, double blMulti); 
  virtual ~ParsimonySPR(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune() ;

  virtual void setOwningChain(Chain *_chain) {chain = _chain;}  
  AbstractProposal* clone() const; 
  
protected: 
  Chain *chain; 
  double parsWarp; 
  double blMulti;   
  Path path; 
}; 
