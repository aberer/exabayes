/**
   @file PriorBelief.hpp

   @brief top level object that represents the prior probability of a
   chain
 */ 


#ifndef _PRIORMANAGER_H
#define _PRIORMANAGER_H

#include <cassert>
#include <memory>
#include <vector>
#include <iostream>

#include "TreeAln.hpp"
#include "Priors.hpp"
#include "GlobalVariables.hpp"
#include "RandomVariable.hpp"


class PriorBelief
{
public:
  PriorBelief(const TreeAln &traln, const vector<RandomVariable>  &variables);

  void accept()  { lnPrior += lnPriorRatio;  lnPriorRatio = 0; }  
  void reject() {lnPriorRatio = 0; }
  double getLnPrior () const {return lnPrior; } 
  double getLnPriorRatio() const {return lnPriorRatio; }
  void addToRatio(double val)  { lnPriorRatio += val; }
  void accountForFracChange(const TreeAln &traln, int model, const vector<double> &oldFc, const vector<double> &newFcs, double lambda )  ; 
  void updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, shared_ptr<AbstractPrior> brPr) ; 
  void verifyPrior(const TreeAln &traln) const ;  
  void reinitPrior(const TreeAln &traln ) {lnPrior = scoreEverything(traln); lnPriorRatio = 0; }
  shared_ptr<AbstractPrior> getBranchLengthPrior() const; 

private: 
  double scoreEverything(const TreeAln &traln) const ; 

  double lnPrior; 
  double lnPriorRatio; 

  // bad, but it's the only possibility
  const vector<bool> needsFranchChangeAccounting; // does changing the fracchange influence the prior of this partition? 
  const vector<double> lambdaOfPartition; 
  
  const vector<RandomVariable> variables; 

}; 

#endif
