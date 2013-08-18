#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include <string>
#include <vector>


#include "Category.hpp"
#include "axml.h"
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
  AbstractProposal()		
  {
    needsFullTraversal = true; 
    inSetExecution = false; 
  } 
  AbstractProposal( const AbstractProposal& rhs); 
  AbstractProposal& operator=(const AbstractProposal &rhs) = delete;  
  virtual ~AbstractProposal(){}

  // you MUST implement all virtual methods in your derived
  // proposal. Here, the signatures are set to 0, this must not
  // be the case in the derived proposal. This 0 keyword makes it
  // impossible to create an instance of AbstractProposal and forces
  // you to implement these methods, when you derive from
  // AbstractProposal.

  /**
     @brief determines the proposal, applies it to the tree / model, updates prior and hastings ratio
   */ 
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) = 0; 
  /**
     @brief evaluates the proposal 
     @todo remove the prior, we should not need it here 
   */ 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln) = 0; 
  /** 
      @brief resets the tree to its previous state; corrects the prior, if necessary (@todo is this the case?)
   */ 
  virtual void resetState(TreeAln &traln) = 0 ; 
  /** 
      @brief tunes proposal parameters, if available  
   */ 
  virtual void autotune() = 0  ;

  virtual AbstractProposal* clone() const = 0;  
  /** 
      @brief gets the relative weight of this proposal 
   */ 
  double getRelativeWeight() const { return relativeWeight; }
  /**
     @brief sets the relative weight of this proposal 
   */ 
  void setRelativeWeight(double tmp) { relativeWeight = tmp; }
  /** 
      @brief gets the category 
   */ 
  Category getCategory() const {return category; }
  /** 
      @brief gets the name of the proposal 
   */ 
  std::string getName() const {return name; }
  /** 
      @brief inform proposal about acceptance
   */ 
  void accept() {sctr.accept();}
  /**
     @brief inform proposal about rejection 
   */ 
  void reject() {sctr.reject();}
  /** 
      @brief gets the success counter 
   */ 
  const SuccessCounter& getSCtr()  const { return sctr; }
  /** 
      @brief gets the number of proposal invocation, since it has been tune the last time 
   */ 
  int  getNumCallSinceTuning() const { return sctr.getRecentlySeen(); }
  /** 
      @brief add a parameter to be integrated over to the proposal 
   */ 
  void addPrimaryParameter(std::unique_ptr<AbstractParameter> var) {primaryParameters.push_back(std::move(var)) ; }
  /** 
      @brief add a parameter  that is integrated over as a by-product of this proposal 
   */ 
  void addSecondaryParameter(std::unique_ptr<AbstractParameter> var) {secondaryParameters.push_back(std::move(var)) ; }
  /** 
      @brief indicates whether this proposal needs a full traversal 
   */ 
  bool isNeedsFullTraversal() const {return needsFullTraversal; }
  std::vector<AbstractParameter*> getBranchLengthsParameterView() const ; 

  /** 
      @brief update the hatsings
      @param valToAdd is the absolute proposal ratio that shall be added 
   */ 
  static void updateHastings(double &hastings, double valToAdd, std::string whoDoneIt); 
  

  std::ostream& printNamePartitions(std::ostream &out); 
  std::ostream& printShort(std::ostream &out)  const ;

  std::vector<AbstractParameter*> getPrimaryParameterView() const; 
  std::vector<AbstractParameter*> getSecondaryParameterView() const ; 
  /** 
      @brief prepare for set execution (only relevant for branch length + node slider)
   */ 
  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  = 0; 

  // we need to implement these 
  void writeToCheckpoint( std::ostream &out)  const; 
  void readFromCheckpoint( std::istream &in ); 
  
  void setPreparedBranch(Branch b ) {preparedBranch = b;  }
  // this is very bad =( 
  void setOtherPreparedBranch(Branch b){preparedOtherBranch = b; }

  void setInSetExecution(bool exec) { inSetExecution = exec;  }
  /** 
      @brief writes proposal specific (tuned) parameters
   */ 
  virtual void writeToCheckpointCore(std::ostream &out)const  = 0 ;  
  /** 
      @brief reads proposal specific (tuned) parameters 
   */ 
  virtual void readFromCheckpointCore(std::istream &in) = 0; 

protected:   
  std::string name;   
  SuccessCounter sctr; 
  Category category; 
  std::vector<std::unique_ptr<AbstractParameter> > primaryParameters; // it is the  primary purpose of this proposal to integrate over these parameters (in most cases only 1) 
  std::vector<std::unique_ptr<AbstractParameter> > secondaryParameters;  // as a by-product also these random variables are changed 
  double relativeWeight; 
  bool needsFullTraversal; 
  bool inSetExecution;

  // meh 
  Branch preparedBranch; 
  Branch preparedOtherBranch; 

  friend std::ostream&  operator<< ( std::ostream& out , AbstractProposal* rhs); 
}; 

#endif
