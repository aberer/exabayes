/**
   @file LnlRestorer.hpp 

   @brief allows to restore the lnl array state before evaluating a proposal 
   
   Note that this doubles the overall memory consumption
 */

#ifndef _LNL_RESTORER_H
#define   _LNL_RESTORER_H

#include "axml.h"
#include "TreeAln.hpp"

typedef node* nodep;  
class Chain;

class LnlRestorer
{
public: 
  /** @brief initializes the arrays (expensive memory-wise) */
  LnlRestorer(Chain *chain);

  ~LnlRestorer();

  /**@brief  resets the restorer, s.t. it is consistent with the current tree (and can restore it later) */ 
  void resetRestorer(TreeAln &traln); 

  /// @brief restores the original Chain of the tree 
  void restoreArrays(TreeAln &traln);  

  /** @brief save any arrays that are recomputed, if they have not   already been flipped */ 
  void traverseAndSwitchIfNecessary(TreeAln &traln, nodep virtualRoot, int model, bool fullTraversal); 

private:   
  void storeOrientation(TreeAln &traln); 
  void loadOrientation(TreeAln &traln);
  void swapArray(TreeAln& traln, int number, int model); 

  int numPart; 
  int numTax; 
  
  // the original ptrs to the allocated stuff; only for cleanup
  // double ***reserveArraysPtr; 

  vector<double> partitionLnl; 

  int modelEvaluated; 
  double ***reserveArrays; 
  int* orientation; 
  bool *wasSwitched; 
  nat **partitionScaler;   
  double prevLnl; 		// for DEBUG
}; 





#endif
