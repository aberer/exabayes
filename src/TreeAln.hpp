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
  explicit TreeAln();
  ~TreeAln();
  TreeAln& operator=( TreeAln &rhs); 
  
  void initRevMat(int model); 
  void unlinkTree();
  nodeptr getUnhookedNode(int number);
  void discretizeGamma(int model); 

  void initializeFromByteFile(char *bytefile); 
  

  /**
     @brief gets the tree length (internal representation)
  */   
  double getTreeLength(); 

  // save setters 
  double setFrequencySave(double newValue, int model, int position ); 
  double setSubstSave(double newValue, int model, int position); 
  double setBranchLengthSave(double newValue, int model, nodeptr p); 
  double setAlphaSave(double newValue, int model); 
  
  void setTr(tree *newTr){ tr = newTr; }
#if HAVE_PLL != 0 
  void setPartitionList(partitionList *pl) { partitions = pl; }
#endif
  

  // getters 
  pInfo* getPartition(int model);
  tree* getTr(){return tr;}
  int getNumBranches(); 
  int getNumberOfPartitions();   
  boolean& accessExecModel(int model); 
  double& accessPartitionLH(int model); 

  // getters for various parameters. We should use them to reduce the
  // direct acceses to for instance branch lengths. Makes our code
  // saver and if we want to do something globally (e.g., print), then
  // it can simply be added to these functions.
  double getBranchLength(nodeptr p, int model){return p->z[model] ;  }
  double getSubstRate(int model, int position){pInfo *partition = getPartition(model); return partition->substRates[position]; }
  double getFrequency(int model, int position){pInfo *partition = getPartition(model); return partition->frequencies[position]; }
  double getAlpha(int model){pInfo *partition = getPartition(model) ; return partition->alpha; } 


  bool isTipNode(nodeptr p){return isTip(p->number, getTr()->mxtips );}

#if HAVE_PLL != 0
  partitionList* getPartitionsPtr(){ return partitions; } 
#endif


  // those are the constraints for various parameters 
  static const double zMin, zMax, 
    rateMax, rateMin, 
    alphaMin, alphaMax,
    freqMin; 
  
  static const double initBL;  	// init values 

  friend ostream& operator<< (ostream& out,  TreeAln&  traln);


private: 
  double getTreeLengthHelper(nodeptr p);
  
  void initDefault();
  tree* tr;		// TODO replace with an object for cleanup 

#if HAVE_PLL != 0
  // horrible hacks, that we cannot get rid of before  upgrading to more recent versions of the PLL 
  partitionList* partitions; 
  void initializeTreePLL();
  void initializePartitionsPLL(char *bytefile, double ***empFreq);
#endif  

};  


#endif
