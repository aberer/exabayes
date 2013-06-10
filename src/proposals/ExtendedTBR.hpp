#include "axml.h"

#include "Path.hpp"
#include "AbstractProposal.hpp"

class ExtendedTBR : public AbstractProposal
{
public: 
  ExtendedTBR( double _extensionProb, double _multiplier); 
  virtual ~ExtendedTBR()  { }

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand); 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior); 
  virtual void resetState(TreeAln& traln, PriorBelief &prior); 
  virtual void autotune();

  virtual AbstractProposal* clone() const; 

private: 
  void drawPaths(TreeAln &traln, Randomness &rand); 
  void executeTBR(TreeAln & traln); 

  double extensionProbability; 
  
  Path modifiedPath1; 
  Path modifiedPath2; 
  double multiplier; 
}; 
