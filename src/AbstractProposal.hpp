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

  double getRelativeWeight() const { return relativeWeight; }
  void setRelativeWeight(double tmp) { relativeWeight = tmp; }
  
  category_t getCategory() const {return category; }
  string getName() const {return name; }
  
  void accept() {sctr.accept();}
  void reject() {sctr.reject();}
  
  const SuccessCounter& getSCtr()  const { return sctr; }
  
  int  getNumCallSinceTuning() const { return sctr.getRecentlySeen(); }

  void addPrimVar(RandomVariable var) {primVar.push_back(var) ; }
  void addSecVar(RandomVariable var) {secVar.push_back(var) ; }

  friend ostream&  operator<< ( ostream& out , const unique_ptr<AbstractProposal> &rhs)
  {
    out << rhs->name <<  " primarily modifying " ; 
    for(auto r : rhs->primVar)
      out << r << ",\t"  ; 

    if(not rhs->secVar.empty() )
      {
	out << "\tand also modifying " ; 
	for(auto r : rhs->secVar ) 
	  out << r << ",\t" ; 
      }

    return out; 
  }

  ostream& printNamePartitions(ostream &out)
  {
    out << name  << "(" ; 
    assert(primVar.size() == 1); 
    bool isFirst= true; 
    for (auto v : primVar[0].getPartitions()) 
      {
	if( not isFirst)
	  out << ","; 
	else 
	  isFirst = false; 
	out << v ; 
      }
    out << ")" ; 
    return out; 
  }

  ostream& printShort(ostream &out) 
  {
    out << name << "( " ;  
    
    bool isFirst = true; 
    for(auto &v : primVar)
      {
	if(not isFirst)
	  out << ","; 
	else 
	  isFirst = false; 
	v.printShort(out); 
      }

    if(secVar.size() > 0)
      {
	out << ";"; 
	isFirst = true; 
	for(auto &v : secVar)
	  {
	    if(not isFirst)
	      out << ","; 
	    else 
	      isFirst = false; 
	    v.printShort(out); 
	  }
      }
    out << " )"; 
    return out; 
  }

protected:   
  string name;   
  SuccessCounter sctr; 
  category_t category; 
  
  vector<RandomVariable> primVar; // it is the  primary purpose of this proposal to integrate over these parameters (in most cases only 1) 
  vector<RandomVariable> secVar;  // as a by-product also these random variables are changed 

  double relativeWeight; 
}; 

#endif
