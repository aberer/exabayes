/** 
    @file ParsimonySPR.hpp
    
    @brief implements a parsimony-biased SPR move similar to MrBayes. 

    @notice: MrBayes also reweights site patters -- currently we do
    not do that.
 */ 

#include "axml.h"
#include "AbstractProposal.hpp"
#include "Path.hpp"
#include "SprMove.hpp"

class ParsimonySPR : public AbstractProposal
{
public: 
  ParsimonySPR(double relativeWeight, double parsWarp, double blMulti); 
  virtual ~ParsimonySPR(){}

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) ; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) ; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) ; 
  virtual void autotune() ;

  AbstractProposal* clone() const; 
  
protected: 
  double parsWarp; 
  double blMulti;   
  Path path; 

  void determineSprPath(TreeAln& traln, Randomness &rand, double &hastings, PriorBelief &prior ); 
  
  SprMove move; 
  
}; 
