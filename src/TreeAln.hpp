#pragma once 


#include "axml.h"
#include "bayes.h"



/** 
    @brief mostly wraps the legacy tr and partition list 
 */

class TreeAln
{
public: 
  TreeAln(tree *tr); 
  ~TreeAln();
  TreeAln& operator=( TreeAln &rhs); 
  TreeAln(const TreeAln &rhs); 

  void initRevMat(int model); 
  void unlinkTree();
  nodeptr getUnhookedNode(int number);

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
  void initFromTree(tree *tr); 

  tree * tr;
#if HAVE_PLL == 1   
  partitionList *partitions; 
#endif  

};  

