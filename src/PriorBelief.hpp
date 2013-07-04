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
#include "parameters/AbstractParameter.hpp"

class PriorBelief
{
public:
  PriorBelief();
  
  void initialize(const TreeAln &traln, std::vector<AbstractParameter*> variables); 
  
  void accept()  { assert(wasInitialized); lnPrior += lnPriorRatio;  lnPriorRatio = 0; }  
  void reject() { assert(wasInitialized) ; lnPriorRatio = 0; }
  double getLnPrior () const {assert(wasInitialized); return lnPrior; } 
  double getLnPriorRatio() const {assert(wasInitialized) ; return lnPriorRatio; }
  void addToRatio(double val)  { assert(wasInitialized) ;  lnPriorRatio += val; }

  void accountForFracChange(const TreeAln &traln, const std::vector<double> &oldFc, const std::vector<double> &newFcs,  const std::vector<AbstractPrior*> &blPriors)  ; 

  void updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, AbstractPrior* brPr) ; 

  void verifyPrior(const TreeAln &traln, std::vector<AbstractParameter*> variables) const ;  

private: 
  double scoreEverything(const TreeAln &traln, std::vector<AbstractParameter*> variables) const ; 
  
  // having an internal state actually defies the logic of the randomVariables being external 
  double lnPrior; 
  double lnPriorRatio; 
  bool wasInitialized; 
}; 



#endif
