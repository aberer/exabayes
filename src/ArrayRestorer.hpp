/**
   @file ArrayRestorer.hpp 

   @brief allows to restore the lnl array state before evaluating a proposal 
   
   Note that this doubles the overall memory consumption
 */

#ifndef _LNL_RESTORER_H
#define   _LNL_RESTORER_H

#include "axml.h"
#include "TreeAln.hpp"

#include <unordered_set>


class ArrayRestorer
{
public: 
  ArrayRestorer(const TreeAln& traln); 
  ~ArrayRestorer();

  /** 
      @brief resets all information about switched arrays
   */ 
  void resetRestorer(const TreeAln &traln); 
  /** 
      @brief restore arrays for the given list of partitions  
   */ 
  void restoreSomePartitions(TreeAln &traln, std::vector<nat> partitions); 
  /** 
      @brief restore all arrays that have been switched
   */ 
  void restoreArrays(TreeAln &traln);  
  /** 
      @brief traverses a subtree and switches arrays, where necessary 
      @param virtualRoot the root of the subtree 
   */ 
  void traverseAndSwitchIfNecessary(TreeAln &traln, nodeptr virtualRoot, std::vector<nat> model, bool fullTraversal); 
  /** 
      @brief traverses the entire tree and switches arrays, where necessary 
   */ 
  void toplevelSwitch(TreeAln &traln, Branch virtualRoot, std::vector<nat> models, bool fullTraversal); 

private:   
  ArrayRestorer(ArrayRestorer &rhs)= delete ; 
  ArrayRestorer& operator=(ArrayRestorer &rhs)= delete ; 

  void storeOrientation(const TreeAln &traln); 
  void loadOrientation(TreeAln &traln);
  void swapArray(TreeAln& traln, int number, std::vector<nat> model); 

  nat numPart; 
  nat numTax; 

  std::unordered_set<nat> modelsEvaluated; 
  std::vector<std::vector<double*> > reservePerPartition; 
  std::vector<std::vector<bool>> wasSwitchedByPartition; /// indicates for each partition, if an array was switched 
  std::vector<nat> orientation;
  std::vector<std::vector<nat> > partitionScaler;   
}; 

#endif
