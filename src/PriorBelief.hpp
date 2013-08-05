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
#include "priors/AbstractPrior.hpp"
#include "GlobalVariables.hpp"
#include "parameters/AbstractParameter.hpp"

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
  void accountForFracChange(const TreeAln &traln, const std::vector<double> &oldFc, const std::vector<double> &newFcs, 
					 const std::vector<AbstractParameter*> &affectedBlParams )  ; 
  /** 
      @brief updates the branch length prior 
      @notice the reason, we have a specific function for that is to avoid some conversions back and forth with the internal representation 
      @param oldInternalZ the old branch length in internal representation 
      @param newInternalZ the new branch length in internal representation 
   */ 
  void updateBranchLengthPrior(const TreeAln &traln , double oldInternalZ,double newInternalZ, const AbstractParameter *param) ; 
  /** 
      @brief verifies the prior 
   */ 
  void verifyPrior(const TreeAln &traln, std::vector<AbstractParameter*> variables) const ;  

  ///////////////
  // OBSERVERS //
  ///////////////
  double getLnPrior () const {assert(wasInitialized); return lnPrior; } 
  double getLnPriorRatio() const {assert(wasInitialized) ; return lnPriorRatio; }

private: 
  double scoreEverything(const TreeAln &traln, const std::vector<AbstractParameter*> &variables) const ; 
  
  // having an internal state actually defies the logic of the randomVariables being external 
  double lnPrior; 
  double lnPriorRatio; 
  bool wasInitialized; 
}; 



#endif
