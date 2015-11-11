/**
   @file PriorBelief.hpp

   @brief top level object that represents the prior probability of a
   chain
 */ 


#ifndef _PRIORMANAGER_H
#define _PRIORMANAGER_H

#include <cassert>
#include <cmath>
#include <memory>
#include <vector>
#include <iostream>

#include "GlobalVariables.hpp"

class AbstractPrior; 
class AbstractParameter; 
class TreeAln;
class PriorBelief
{
public:
  PriorBelief();  
  /** 
      @brief intializes and scores the prior each parameter    
   */ 
  void initialize(const TreeAln &traln, const std::vector<AbstractParameter*> &variables);   
  /** 
      @brief informs the prior about acceptance of the new state    
   */ 
  void accept()  { assert(wasInitialized); lnPrior += lnPriorRatio;  lnPriorRatio = 0; }  
  /** 
      @brief informs the prior about rejection of the new state 
   */
  void reject() { assert(wasInitialized) ; lnPriorRatio = 0; }
  /** 
      @brief adds a (logarithmic!) value to the prior ratio
   */ 
  void addToRatio(double val)  { assert(wasInitialized) ;  lnPriorRatio += val; }
  /** 
      @brief accounts for branch length prior changes due to either
      substitution rate or state frequencies updates

      @param oldFc old fracchanges for each involved branch length parameter
      @param newFc new fracchanges (after proposal) for each involved branch length parameter 

      @notice indices in the fracchange arrays correspond to
      idOfMyKind of branch length parameters

      @notice this is very pedastrian, but I do not see how to avoid this
   */ 
  void accountForFracChange( TreeAln &traln, const std::vector<double> &oldFc, const std::vector<double> &newFcs, 
			    const std::vector<AbstractParameter*> &affectedBlParams )  ; 
  /** 
      @brief verifies the prior 
   */ 
  void verifyPrior(const TreeAln &traln, std::vector<AbstractParameter*> variables) const ;  

  ///////////////
  // OBSERVERS //
  ///////////////
  double getLnPrior () const {assert(not std::isinf(lnPrior)) ; assert(wasInitialized); return lnPrior; } 
  double getLnPriorRatio() const {assert(wasInitialized) ; return lnPriorRatio; }

private: 
  double scoreEverything(const TreeAln &traln, const std::vector<AbstractParameter*> &variables) const ; 
  
  // having an internal state actually defies the logic of the randomVariables being external 
  double lnPrior; 
  double lnPriorRatio; 
  bool wasInitialized; 
}; 



#endif
