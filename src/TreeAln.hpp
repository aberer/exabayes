#pragma once 


#include "axml.h"
#include "bayes.h"



/** @brief
    mostly wraps the legacy tr
 */

class TreeAln
{
public: 
  TreeAln(tree *tr ,bool allocateStructures); 
  TreeAln(const TreeAln &rhs); 
  ~TreeAln();
  TreeAln& operator=( TreeAln &rhs); 
  


  void initRevMat(int model); 


  // getters 
  pInfo* getPartition(int model);
  tree* getTr(){return tr;}
  int getNumBranches(); 
  int getNumberOfPartitions();   
  int& accessExecModel(int model); 
  double& accessPartitionLH(int model); 

#if HAVE_PLL == 1 
  partitionList* getPartitionsPtr(){return partitions; } 
#endif

private: 
  void initDefault();
  void initFromTree(tree *tr) ; 
  void unlinkTree();
  nodeptr getUnhookedNode(int number);

  tree * tr;
#if HAVE_PLL == 1   
  partitionList *partitions; 
#endif  

};  



