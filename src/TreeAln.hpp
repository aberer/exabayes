#ifndef _TREEALN_H
#define _TREEALN_H

#include <vector>
#include <iostream>
#include <memory>


#include "Priors.hpp"
#include "Branch.hpp"
#include "axml.h"
#include "Randomness.hpp"

#include "TreePrinter.hpp"

nat numStateToNumInTriangleMatrix(int numStates) ; 

/** 
    @brief mostly wraps the legacy tr and partition list     
*/

class TreeAln
{

public: 
   /////////////////
   // life cycle  //
   /////////////////
  TreeAln();
  ~TreeAln();
  TreeAln& operator=( TreeAln &rhs); 
  TreeAln(const TreeAln &tmp) = delete ; 


  /////////////////////////////////////
  //           OBSERVERS             //
  /////////////////////////////////////
  pInfo* getPartition(int model) const;
  tree* getTr() {return &tr;}
  const tree* getTr() const{return &tr; }
  int getNumBranches() const; 
  int getNumberOfPartitions() const;   
  nat getNumberOfTaxa() const {return getTr()->mxtips; }
  nat getNumberOfNodes() const { nat numTax = getNumberOfTaxa(); return 2 * numTax - 3 ;  } // excluding the virtual root 

  

  ///////////////
  // MODIFIERS //
  ///////////////
  void setFrequencies(const std::vector<double> &values, int model);
  void setRevMat(const std::vector<double> &values, int model);
  void setAlpha(double alpha,  int model);   
  void setBranch(const Branch& b);   

   

  /** @notice returns the revmat, s.t. rates sum up to 1 */ 
  std::vector<double> getRevMat(int model) const ;   
  std::vector<double> getFrequencies(int model) const; 
  double getAlpha(int model) const {pInfo *partition = getPartition(model) ; return partition->alpha; } 
  bool isTipNode(nodeptr p) const {return isTip(p->number, getTr()->mxtips );}
  Branch getBranch(nodeptr p) const { return Branch(p->number, p->back->number, p->z[0]); }
  nodeptr getNode(nat elem) const ; 
  std::vector<Branch> extractBranches() const ; 
  double getTreeLengthExpensive() const;

#if HAVE_PLL != 0 
  void setPartitionList(partitionList *pl) { partitions = *pl; }
  partitionList* getPartitionsPtr()  { return &partitions; }   
#endif  


  // MODIFIERS 
  void clipNode(nodeptr p, nodeptr q, double z);   
  void clipNodeDefault(nodeptr p, nodeptr q); 
  void enableParsimony();
  void unlinkTree();
  void initializeFromByteFile(std::string  byteFileName);   

  // TODO private 
  void initRevMat(int model); 
  void discretizeGamma(int model); 

  nodeptr getUnhookedNode(int number);

  // TODO 
  // void setFrequenciesBounded(std::vector<double> &newValues, int model ); 
  // void setRevMatBounded(std::vector<doub2le> &newValues, int model); 
  // void setBranchLengthBounded(double &newValue, int model, nodeptr p); 
  // void setAlphaBounded(double &newValue, int model); 


  ///////////////
  // observers //
  ///////////////
  boolean& accessExecModel(int model); 
  double& accessPartitionLH(int model); 
  int accessExecModel(int model) const; 
  double accessPartitionLH(int model) const ; 

  static const double zZero;   
  static const double initBL;  	// init values 


private: 			// METHODS
#if HAVE_PLL != 0
  void initializeTreePLL(std::string byteFileName);
  void initializePartitionsPLL(std::string byteFileName, double ***empFreq, bool multiBranch);
#endif  
  double getTreeLengthHelper(nodeptr p) const;
  void extractHelper( nodeptr p , std::vector<Branch> &result, bool isStart) const ; 
  void initDefault();



private: 			// ATTRIBUTES 
#if HAVE_PLL != 0 
  // horrible hacks, that we cannot get rid of before  upgrading to more recent versions of the PLL 
  partitionList partitions; 
#endif
  tree tr;		// TODO replace with an object for cleanup   
  bool parsimonyEnabled;   



  // friends 
  friend std::ostream& operator<< (std::ostream& out,  TreeAln&  traln);

  //////////////////
  // EXPERIMENTAL //
  //////////////////
public: 
  // BEGIN collapse test 
  void collapseBranch(Branch b); 
  bool isCollapsed(Branch b) ; 
  void setBranchLengthUnsafe(Branch b ) ; 
  // END

  bool revMatIsImmutable(int model) const ; 

  bool usingPerSiteRates()  { return tr.rateHetModel == GAMMA;  } 
  void enablePerSiteRates() {  tr.rateHetModel = CAT; } 

};  


#endif
