#ifndef _TREEALN_H
#define _TREEALN_H

#include <vector>
#include <iostream>
#include <memory>

using namespace std; 

#include "axml.h"

// #include "LnlRestorer.hpp"
class LnlRestorer; 

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

  void initializeFromByteFile(string  byteFileName); 
  

  void enableParsimony();

  /**
     @brief gets the tree length (internal representation)
  */   
  double getTreeLength() const ; 

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
  pInfo* getPartition(int model) const;
  tree* getTr() const {return tr;}
  int getNumBranches() const; 
  int getNumberOfPartitions() const;   
  boolean& accessExecModel(int model); 
  double& accessPartitionLH(int model); 

  // getters for various parameters. We should use them to reduce the
  // direct acceses to for instance branch lengths. Makes our code
  // saver and if we want to do something globally (e.g., print), then
  // it can simply be added to these functions.
  double getBranchLength(nodeptr p, int model) const {return p->z[model] ;  }
  double getSubstRate(int model, int position) const { return getPartition(model)->substRates[position]; }
  double getFrequency(int model, int position) const {pInfo *partition = getPartition(model); return partition->frequencies[position]; }
  double getAlpha(int model) const {pInfo *partition = getPartition(model) ; return partition->alpha; } 


  bool isTipNode(nodeptr p) const {return isTip(p->number, getTr()->mxtips );}

#if HAVE_PLL != 0
  partitionList* getPartitionsPtr() const { return partitions; } 
#endif


  // those are the constraints for various parameters 
  static const double zMin, zMax, 
    rateMax, rateMin, 
    alphaMin, alphaMax,
    freqMin; 
  
  static const double initBL;  	// init values 

  friend ostream& operator<< (ostream& out,  TreeAln&  traln);

  void setRestorer(shared_ptr<LnlRestorer> rest){ restorer = rest; }
  auto getRestorer() const ->  shared_ptr<LnlRestorer> {return restorer;  }  


private:   
  double getTreeLengthHelper(nodeptr p) const;
  
  void initDefault();
  tree* tr;		// TODO replace with an object for cleanup 
  
  bool parsimonyEnabled; 
  
  bool branchLengthsFixed; 
  
  shared_ptr<LnlRestorer> restorer; 

#if HAVE_PLL != 0
  // horrible hacks, that we cannot get rid of before  upgrading to more recent versions of the PLL 
  partitionList* partitions; 
  void initializeTreePLL(string byteFileName);
  // void initializeTreePLL();
  void initializePartitionsPLL(string byteFileName, double ***empFreq, bool multiBranch);
#endif  

};  


#endif
