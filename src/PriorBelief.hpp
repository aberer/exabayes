/**
   @file PriorBelief.hpp

   @brief top level object that represents the prior probability of a
   chain
 */ 


#ifndef _PRIORMANAGER_H
#define _PRIORMANAGER_H

using namespace std; 

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

  void setBranchLengthPrior(shared_ptr<AbstractPrior> ptr)  { brPr = ptr ; }
  void setRevMatPrior(shared_ptr<AbstractPrior> ptr)  {revMatPr = ptr; }
  void setRateHetPrior(shared_ptr<AbstractPrior> ptr)  {rateHetPr = ptr; }
  void setStateFreqPrior(shared_ptr<AbstractPrior> ptr)  {stateFreqPr = ptr; }

  double getLogProb(){return logProb;}

  void verify(const TreeAln& traln);

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
