#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include "axml.h"
#include <string>
#include "SuccessCtr.hpp" 
#include "proposalType.h"
#include "Randomness.hpp"
#include "categoryType.h" 
#include "proposalType.h"
#include "PriorManager.hpp"
#include "GlobalVariables.hpp"

class Chain; 

class AbstractProposal
{
public: 
  AbstractProposal(){}
  virtual ~AbstractProposal(){}


  // you MUST implement all virtual methods in your derived
  // proposal. Here, the signatures are set to 0, this must not
  // be the case in the derived proposal. This 0 keyword makes it
  // impossible to create an instance of AbstractProposal and forces
  // you to implement these methods, when you derive from
  // AbstractProposal.
  virtual void applyToState(TreeAln &traln, PriorManager &prior, double &hastings, Randomness &rand) = 0; 
  virtual void evaluateProposal(TreeAln &traln, PriorManager &prior) = 0; 
  virtual void resetState(TreeAln &traln, PriorManager &prior) = 0 ; 
  virtual void autotune() = 0  ;

  // HACK! 
  virtual void setOwningChain( Chain *chain) =  0 ; 

  // getters and setters 
  double getRelativeProbability() const {return relativeProbability; }
  void setRelativeProbability(double prob ){ relativeProbability = prob; }
  
  category_t getCategory() const {return category; }
  proposal_type getPtype() const {return ptype; }
  
  string getName() const {return name; }
  
  void accept() {sctr.accept();}
  void reject() {sctr.reject();}
  
  SuccessCtr getSCtr() const { return sctr; }
  
  int  getNumCallSinceTuning(){ return sctr.getRecentlySeen(); }


protected: 
  string name;   
  SuccessCtr sctr; 
  category_t category; 
  double relativeProbability ; 	// probability relative to category
  proposal_type ptype;  		// TODO will not be needed later  


  // various parameters and stuff 

}; 


#endif
