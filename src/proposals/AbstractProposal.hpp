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
#include "LikelihoodEvaluator.hpp"
#include "TreeRandomizer.hpp"
#include "Checkpointable.hpp"


class AbstractProposal : public Checkpointable
{
public: 
  // for copying the non-trivial types 
  AbstractProposal(){} 
  AbstractProposal( const AbstractProposal& rhs); 
  
  virtual ~AbstractProposal(){}

  // you MUST implement all virtual methods in your derived
  // proposal. Here, the signatures are set to 0, this must not
  // be the case in the derived proposal. This 0 keyword makes it
  // impossible to create an instance of AbstractProposal and forces
  // you to implement these methods, when you derive from
  // AbstractProposal.
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) = 0; 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, PriorBelief &prior) = 0; 
  virtual void resetState(TreeAln &traln, PriorBelief &prior) = 0 ; 
  virtual void autotune() = 0  ;

  virtual AbstractProposal* clone() const = 0;  
  double getRelativeWeight() const { return relativeWeight; }
  void setRelativeWeight(double tmp) { relativeWeight = tmp; }
  Category getCategory() const {return category; }
  std::string getName() const {return name; }
  void accept() {sctr.accept();}
  void reject() {sctr.reject();}
  const SuccessCounter& getSCtr()  const { return sctr; }
  int  getNumCallSinceTuning() const { return sctr.getRecentlySeen(); }
  void addPrimVar(std::unique_ptr<AbstractParameter> var) {primVar.push_back(std::move(var)) ; }
  void addSecVar(std::unique_ptr<AbstractParameter> var) {secVar.push_back(std::move(var)) ; }
  static void updateHastings(double &hastings, double valToAdd, std::string whoDoneIt); 
  friend std::ostream&  operator<< ( std::ostream& out , const std::unique_ptr<AbstractProposal> &rhs); 
  std::ostream& printNamePartitions(std::ostream &out); 
  std::ostream& printShort(std::ostream &out)  const ;

  std::vector<AbstractParameter*> getPrimVar() const; 
  std::vector<AbstractParameter*> getSecVar() const ; 

  void writeToCheckpoint( std::ofstream &out) const ; 
  void readFromCheckpoint( std::ifstream &in ); 

  virtual void writeToCheckpointCore(std::ofstream &out) const = 0 ;  
  virtual void readFromCheckpointCore(std::ifstream &in) = 0; 

protected:   
  std::string name;   
  SuccessCounter sctr; 
  Category category; 

  // will be a unique_ptr later 
  std::vector<std::unique_ptr<AbstractParameter> > primVar; // it is the  primary purpose of this proposal to integrate over these parameters (in most cases only 1) 
  std::vector<std::unique_ptr<AbstractParameter> > secVar;  // as a by-product also these random variables are changed 

  double relativeWeight; 
}; 

#endif
