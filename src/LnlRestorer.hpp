/**
   @file LnlRestorer.hpp 

   @brief allows to restore the lnl array state before evaluating a proposal 
   
   Note that this doubles the overall memory consumption
 */

#ifndef _LNL_RESTORER_H
#define   _LNL_RESTORER_H

#include "axml.h"
#include "bayes.h"

typedef node* nodep;  


class LnlRestorer
{
public: 
  /** @brief initializes the arrays (expensive memory-wise) */
  LnlRestorer(Chain *chain);

  ~LnlRestorer();

  /**@brief  resets the restorer, s.t. it is consistent with the current tree (and can restore it later) */ 
  void resetRestorer(); 

  /// @brief restores the original Chain of the tree 
  void restore();  

  /** @brief save any arrays that are recomputed, if they have not   already been flipped */ 
  void traverseAndSwitchIfNecessary(nodep virtualRoot, int model, bool fullTraversal); 


private:   
  void storeOrientation(); 
  void loadOrientation();
  void swapArray(int number, int model); 
  
  // the original ptrs to the allocated stuff; only for cleanup
  // double ***reserveArraysPtr; 

  int modelEvaluated; 
  Chain *chain; 
  double ***reserveArrays; 
  int* orientation; 
  bool *wasSwitched; 
  nat **partitionScaler;   
  double prevLnl; 		// for DEBUG
}; 





#endif
