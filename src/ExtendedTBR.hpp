#include "bayes.h"
#include "axml.h"

#include "Path.hpp"
#include "AbstractProposal.hpp"

class Chain; 


class ExtendedTBR : public AbstractProposal
{
public: 
  ExtendedTBR( Chain *chain, double relativeProbability, double extensionProb, double multiplier);
  virtual ~ExtendedTBR(){};

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand); 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior); 
  virtual void resetState(TreeAln& traln, PriorManager &prior); 
  virtual void autotune();
  
  virtual void setOwningChain(Chain *chain){}
  
private: 
  void drawPaths(TreeAln &traln, Randomness &rand); 
  void executeTBR(TreeAln & traln); 
  

  Chain *chain; 
  double extensionProbability; 
  
  Path modifiedPath1; 
  Path modifiedPath2; 
  double multiplier; 

}; 
