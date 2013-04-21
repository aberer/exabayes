#ifndef _TREEALN_H
#define _TREEALN_H


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
  void discretizeGamma(int model); 
  

  // save setters 
  double setFrequencySave(double newValue, int model, int position ); 
  double setSubstSave(double newValue, int model, int position); 
  double setBranchLengthSave(double newValue, int model, nodeptr p); 
  double setAlphaSave(double newValue, int model); 

  // getters 
  pInfo* getPartition(int model);
  tree* getTr(){return tr;}
  int getNumBranches(); 
  int getNumberOfPartitions();   
  int& accessExecModel(int model); 
  double& accessPartitionLH(int model); 

  // getters for various parameters. We should use them to reduce the
  // direct acceses to for instance branch lengths. Makes our code
  // saver and if we want to do something globally (e.g., print), then
  // it can simply be added to these functions.
  double getBranchLength(nodeptr p, int model){return p->z[model] ;  }
  double getSubstRate(int model, int position){pInfo *partition = getPartition(model); return partition->substRates[position]; }
  double getFrequency(int model, int position){pInfo *partition = getPartition(model); return partition->frequencies[position]; }
  double getAlpha(int model){pInfo *partition = getPartition(model) ; return partition->alpha; } 


#if HAVE_PLL != 0
  partitionList* getPartitionsPtr(){return partitions; } 
#endif


  // those are the constraints for various parameters 
  static const double zMin, zMax, 
    rateMax, rateMin, 
    alphaMin, alphaMax,
    freqMin; 



private: 
  void initDefault();
  void initFromTree(tree *tr); 

  tree * tr;
#if HAVE_PLL != 0
  partitionList *partitions; 
#endif  

};  


#endif
