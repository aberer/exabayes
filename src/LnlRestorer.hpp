/**
   @file LnlRestorer.hpp 

   @brief allows to restore the lnl array state before evaluating a proposal 
   
   Note that this doubles the overall memory consumption
 */

#ifndef _LNL_RESTORER_H
#define   _LNL_RESTORER_H

#include "axml.h"
#include "TreeAln.hpp"

#define ALL_MODELS  -1  

class LnlRestorer
{
public: 
  /** @brief initializes the arrays (expensive memory-wise) */
  LnlRestorer(TreeAln& traln); 

  ~LnlRestorer();

  /**@brief  resets the restorer, s.t. it is consistent with the current tree (and can restore it later) */ 
  void resetRestorer(const TreeAln &traln); 

  /// @brief restores the original Chain of the tree 
  void restoreArrays(TreeAln &traln);  

  /** @brief save any arrays that are recomputed, if they have not   already been flipped */ 
  void traverseAndSwitchIfNecessary(TreeAln &traln, nodeptr virtualRoot, int model, bool fullTraversal); 

private:   
  LnlRestorer(LnlRestorer &rhs)= delete ; 
  LnlRestorer& operator=(LnlRestorer &rhs)= delete ; 

  void storeOrientation(const TreeAln &traln); 
  void loadOrientation(TreeAln &traln);
  void swapArray(TreeAln& traln, int number, int model); 

  int numPart; 
  int numTax; 

  std::vector<double> partitionLnl; 

  int modelEvaluated; 
  double ***reserveArrays; 
  std::vector<bool> wasSwitched; 
  std::vector<nat> orientation; 
  std::vector<std::vector<nat>> partitionScaler;   
  double prevLnl; 		// for DEBUG

}; 

#endif
