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
#include "GlobalVariables.hpp"


class PriorBelief
{
public: 
  PriorBelief(); 
  ~PriorBelief(){};

  // TODO support topology 

  void updateBranchLengthAllByFactor(double factor)
  {
    cout << "branch lengths prior: "  << lnBrProb << " => " << lnBrProb * factor << endl; 
    lnBrProb *= factor; 
  }


  /** @brief updates the branch lengths prior with a new value (and cancels out the old one) */ 
  void updateBranchLength(double oldValue, double newValue); 

  /** @brief updates the revMat prior with a new value (and cancels out the old one) */ 
  void updateRevMat( vector<double> oldValues, vector<double> newValues); 

  /** @brief updates the state frequency prior with a new value (and cancels out the old one) */ 
  void updateFreq(vector<double> oldValues, vector<double> newValues); 

  /** @brief updates the rate heterogeneity prior with a new value (and cancels out the old one) */ 
  void updateRateHet(double oldValue, double newValue); 

  void initPrior(const TreeAln &traln);

  
  /** @brief sets the bl prior */ 
  void setBranchLengthPrior(shared_ptr<AbstractPrior> ptr)  { assertEmpty(brPr) ;  ; brPr = ptr ; }
  
  /** @brief sets the revmat  prior */ 
  void setRevMatPrior(shared_ptr<AbstractPrior> ptr)  {assertEmpty(revMatPr);  revMatPr = ptr; }

  /** @brief sets the rate het prior */ 
  void setRateHetPrior(shared_ptr<AbstractPrior> ptr)  { assertEmpty(rateHetPr);  rateHetPr = ptr; }
  
  /** @brief sets the state freq prior */ 
  void setStateFreqPrior(shared_ptr<AbstractPrior> ptr)  { assertEmpty(stateFreqPr); stateFreqPr = ptr; }

  /** @brief sets the topology prior */ 
  void setTopologyPrior(shared_ptr<AbstractPrior> ptr) {assertEmpty(topologyPr) ; topologyPr = ptr; }
  
  
  /** @brief gets the logarithmic probabililty of all parameters   */ 
  double getLogProb() const {return lnBrProb + lnRevMatProb + lnRateHetProb + lnStatFreqProb;}

  /** @brief verifies that the currently stored probability is correct */ 
  void verify(const TreeAln& traln);
 
  /**
     @brief adds some standard priors, in case the user has not
     provided specifications for all type of priors
  */ 
  void addStandardPriors(const TreeAln &traln ); 


  shared_ptr<AbstractPrior> getBranchLengthPrior(){ return brPr; }

  // print method 
  friend ostream& operator<<(ostream &out, const PriorBelief& rhs)
  {
    return out
      << "Topology:\t" << rhs.topologyPr.get() << endl 
      << "Branch Lengths:\t" << rhs.brPr.get() << endl  
      << "Reversible Matrix:\t" << rhs.revMatPr.get() << endl 
      << "Rate Heterogeneity:\t" << rhs.rateHetPr.get() << endl 
      << "State Frequencies:\t" << rhs.stateFreqPr.get() << endl; 
  }

  double scoreBranchLengths(const TreeAln &traln); 

  void  rescoreAllBranchLengths(const TreeAln &traln)
  {
#ifdef EFFICIENT
    assert(0); 
#endif
    lnBrProb = scoreBranchLengths(traln); 
  }
  

  bool categoryIsFixed(category_t cat) const
  {
    switch(cat)
      {
      case TOPOLOGY: 
	return believingInFixedTopology(); 
      case BRANCH_LENGTHS: 
	return believingInFixedBranchLengths();
      case SUBSTITUTION_RATES: 
	return believingInFixedRevMat(); 
      case RATE_HETEROGENEITY: 
	return believingInFixedRateHet(); 
      case FREQUENCIES: 
	return believingInFixedStateFreq();
      default: 
	assert(0); 
      }
  }


  vector<double> drawFromPriorByCategory(category_t cat, Randomness &rand)
  {
    if(cat == BRANCH_LENGTHS)
      {
	return brPr->drawFromPrior(rand);
      }
    else 
      {
	assert(0) ; 
	vector<double>tmp; 
	return tmp; 
      }
  }


  // checks if anything is fixed 
  bool believingInFixedBranchLengths() const { return  typeid(*(brPr.get())) == typeid(FixedPrior) ;  } 
  bool believingInFixedRevMat() const { return typeid(*(revMatPr.get())) == typeid(FixedPrior);  }
  bool believingInFixedRateHet() const {return typeid(*(rateHetPr.get())) == typeid(FixedPrior) ; }
  bool believingInFixedStateFreq() const { return typeid(*(stateFreqPr.get())) == typeid(FixedPrior) ;  }
  bool believingInFixedTopology() const {return typeid(*(topologyPr.get())) == typeid(FixedPrior); }

  
private:
  double lnBrProb; 		// prior for branch lengths
  double lnRevMatProb; 		// prior for revmat 
  double lnRateHetProb; 	// prior for rate hetero 
  double lnStatFreqProb; 	// prior for state frequencies

  shared_ptr<AbstractPrior> topologyPr; // this one is still treated special somewhat 
  shared_ptr<AbstractPrior> brPr; 
  shared_ptr<AbstractPrior> revMatPr; 
  shared_ptr<AbstractPrior> rateHetPr; 
  shared_ptr<AbstractPrior> stateFreqPr; 

  double scoreEverything(const TreeAln& traln); 

  double scoreRevMats(const TreeAln &traln); 
  double scoreRateHets(const TreeAln &traln);
  double scoreStateFreqs(const TreeAln &traln); 

  void assertEmpty(shared_ptr<AbstractPrior> prior);
}; 


#endif
