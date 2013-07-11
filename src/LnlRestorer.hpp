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
  LnlRestorer(TreeAln& traln); 
  // ~LnlRestorer();

  void resetRestorer(const TreeAln &traln); 
  void restoreArrays(TreeAln &traln);  
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
  // TODO eventually use unique ptrs there as well 
  std::vector<std::vector<double* > > reservePerPartition; 
  std::vector<bool> wasSwitched; 
  std::vector<nat> orientation; 
  std::vector<std::vector<nat>> partitionScaler;   
  double prevLnl; 		// for DEBUG
}; 

#endif
