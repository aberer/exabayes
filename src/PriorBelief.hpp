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
  PriorBelief(const TreeAln &traln, const vector<RandomVariablePtr> &variables);

  void accept()  { lnPrior += lnPriorRatio;  lnPriorRatio = 0; }  
  void reject() {lnPriorRatio = 0; }
  double getLnPrior () const {return lnPrior; } 
  double getLnPriorRatio() const {return lnPriorRatio; }
  void addToRatio(double val)  { lnPriorRatio += val; }

  void accountForFracChange(const TreeAln &traln, const vector<double> &oldFc, const vector<double> &newFcs,  const vector<shared_ptr<AbstractPrior> > &blPriors)  ; 

  void updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, shared_ptr<AbstractPrior> brPr) ; 
  void verifyPrior(const TreeAln &traln, const vector<RandomVariablePtr> &variables) const ;  
  void reinitPrior(const TreeAln &traln, const vector<RandomVariablePtr> &variables ) {lnPrior = scoreEverything(traln, variables); lnPriorRatio = 0; }

private: 
  double scoreEverything(const TreeAln &traln, const vector<RandomVariablePtr> &variables) const ; 
  
  // having an internal state actually defies the logic of the randomVariables being external 
  double lnPrior; 
  double lnPriorRatio; 
}; 

#endif
