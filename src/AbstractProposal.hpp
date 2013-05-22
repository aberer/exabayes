#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include <string>
#include <vector>

#include "axml.h"
#include "RandomVariable.hpp"
#include "SuccessCounter.hpp" 
#include "Randomness.hpp"
#include "PriorBelief.hpp"
#include "GlobalVariables.hpp"

void updateHastings(double &hastings, double valToAdd, string whoDoneIt); 

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
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) = 0; 
  virtual void evaluateProposal(TreeAln &traln, PriorBelief &prior) = 0; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) = 0 ; 
  virtual void autotune() = 0  ;
  
  virtual void smallPrint(ostream &out) {assert(NOT_IMPLEMENTED); } // TODO 
  virtual void printEverything(ostream &out){assert(NOT_IMPLEMENTED); }  // TODO 

  virtual AbstractProposal* clone() const = 0;  

  virtual double getRelativeWeight() const = 0; 


  category_t getCategory() const {return category; }
  string getName() const {return name; }
  
  void accept() {sctr.accept();}
  void reject() {sctr.reject();}
  
  const SuccessCounter& getSCtr()  const { return sctr; }
  
  int  getNumCallSinceTuning(){ return sctr.getRecentlySeen(); }

  void addRandomVariable(RandomVariable var) {randomVariables.push_back(var) ; }



  friend ostream&  operator<< ( ostream& out , const unique_ptr<AbstractProposal> &rhs)
  {
    out << "proposal " << rhs->name <<  " modifying " ; 
    for(auto r : rhs->randomVariables)
      out << r << ",\t"  ; 
    return out; 
  }



protected: 
  string name;   
  SuccessCounter sctr; 
  category_t category; 
  vector<RandomVariable> randomVariables; // random variables that are integrated by this proposal
}; 

#endif
