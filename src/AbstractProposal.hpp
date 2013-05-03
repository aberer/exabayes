#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include <string>
#include "SuccessCtr.hpp" 
#include "proposalType.h"
#include "Randomness.hpp"
#include "categoryType.h" 
#include "proposalType.h"
#include "PriorManager.hpp"


class Chain; 

class AbstractProposal
{
public: 
  AbstractProposal(){}
  virtual ~AbstractProposal(){}

  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) = 0; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) = 0; 
  virtual void resetState(TreeAln &traln, PriorManager &prior) = 0; 
  virtual void autotune() = 0;

  // HACK! 
  virtual void setOwningChain( Chain *chain) = 0; 
  
  
  // getters and setters 
  double getRelativeProbability(){return relativeProbability; }
  void setRelativeProbability(double prob ){ relativeProbability = prob; }
  
  category_t getCategory(){return category; }
  proposal_type getPtype(){return ptype; }
  
  string getName(){return name; }
  
  void accept() {sctr.accept();}
  void reject() {sctr.reject();}
  
  SuccessCtr getSCtr() const { return sctr; }
  
  bool isTimeToTune(int tuneFreq){return true;  }

protected: 
  string name;   
  SuccessCtr sctr; 
  category_t category; 
  double relativeProbability ; 	// probability relative to category
  proposal_type ptype;  		// TODO will not be needed later  


  // various parameters and stuff 

}; 


#endif
