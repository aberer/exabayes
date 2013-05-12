/**
   @file PriorBelief.hpp

   @brief top level object that represents the prior probability of a
   chain
 */ 


#ifndef _PRIORMANAGER_H
#define _PRIORMANAGER_H

using namespace std; 
#include <cassert>
#include <memory>
#include <vector>
#include <iostream>

#include "TreeAln.hpp"
#include "Priors.hpp"


class PriorBelief
{
public: 
  PriorBelief(); 
  ~PriorBelief(){};

  // TODO support topology 

  void updateBranchLength(double oldValue, double newValue); 
  void updateRevMat( vector<double> oldValues, vector<double> newValues); 
  void updateFreq(vector<double> oldValues, vector<double> newValues); 
  void updateRateHet(double oldValue, double newValue); 

  double scoreEverything(const TreeAln &traln); 
  void initPrior(const TreeAln &traln);

  void setBranchLengthPrior(shared_ptr<AbstractPrior> ptr)  { assertEmpty(brPr) ;  ; brPr = ptr ; }
  void setRevMatPrior(shared_ptr<AbstractPrior> ptr)  {assertEmpty(revMatPr);  revMatPr = ptr; }
  void setRateHetPrior(shared_ptr<AbstractPrior> ptr)  { assertEmpty(rateHetPr);  rateHetPr = ptr; }
  void setStateFreqPrior(shared_ptr<AbstractPrior> ptr)  { assertEmpty(stateFreqPr); stateFreqPr = ptr; }

  double getLogProb() const {return logProb;}

  void verify(const TreeAln& traln);

  void assertEmpty(shared_ptr<AbstractPrior> prior)
  {
    if(prior.get() != nullptr)
      {
	cerr << "Error! Attempted to set prior already specified prior "  << prior.get() << ". Did you specify it twice in the config file?" << endl; 
	assert(0); 
      }
  } 

  /**
     @brief adds some standard priors, in case the user has not provided
     specifications for all type of priors
  */ 
  void addStandardPriors(); 


  friend ostream& operator<<(ostream &out, const PriorBelief& rhs)
  {
    return out << "Branch Lengths:\t" << rhs.brPr.get() << endl  
	       << "Reversible Matrix:\t" << rhs.revMatPr.get() << endl 
	       << "Rate Heterogeneity:\t" << rhs.rateHetPr.get() << endl 
	       << "State Frequencies:\t" << rhs.stateFreqPr.get() << endl; 
  }
  
  
private: 
  double logProb; 

  shared_ptr<AbstractPrior> brPr; 
  shared_ptr<AbstractPrior> revMatPr; 
  shared_ptr<AbstractPrior> rateHetPr; 
  shared_ptr<AbstractPrior> stateFreqPr; 

  double scoreBranchLengths(const TreeAln &traln); 
  double scoreRevMats(const TreeAln &traln); 
  double scoreRateHets(const TreeAln &traln);
  double scoreStateFreqs(const TreeAln &traln); 
}; 


#endif
