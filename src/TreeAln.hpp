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
  /**
     @brief copies the entire state from the rhs to this tree/alignment.
     Important: likelihood arrays are not concerned
  */ 
  TreeAln& operator=( TreeAln &rhs); 
  TreeAln(const TreeAln &tmp) = delete ; 


  /////////////////////////////////////
  //           OBSERVERS             //
  /////////////////////////////////////
  /**
     @brief get the raxml representation of the partition  
   */ 
  pInfo* getPartition(int model) const;
  /** 
      @brief get the internal raxml tree representation 
   */ 
  tree* getTr() {return &tr;}
  /** 
      @brief get a constant internal raxml tree representation  
   */ 
  const tree* getTr() const{return &tr; }
  /** 
      @brief get the number of per-partition branch lengths
   */ 
  int getNumBranches() const; 
  /** 
      @brief get the number of branches in the tree (not counting per-partition branch lengths)
   */ 
  nat getNumberOfBranches() const {return getNumberOfNodes() - 1 ; }
  /** 
      @brief get the number of partitions in the alignment
   */ 
  int getNumberOfPartitions() const;   
  /** 
      @brief get the number of taxa in the tree 
   */ 
  nat getNumberOfTaxa() const {return getTr()->mxtips; }
  /** 
      @brief get the number of nodes in the unrooted tree 
   */ 
  nat getNumberOfNodes() const { nat numTax = getNumberOfTaxa(); return 2 * numTax - 3 ;  } // excluding the virtual root 
  /** 
      @brief get the substitution matrix for partition "model"
   */ 
  std::vector<double> getRevMat(int model) const ;   
  /** 
      @brief gets the state frequencies for partition "model" 
   */ 
  std::vector<double> getFrequencies(int model) const; 
  /**
     @brief gets the alpha parameter of the gamma distribution  
   */ 
  double getAlpha(int model) const {pInfo *partition = getPartition(model) ; return partition->alpha; } 
  /** 
      @brief indicates whether a nodepointer is a tip 
   */ 
  bool isTipNode(nodeptr p) const {return isTip(p->number, getTr()->mxtips );}
  /** 
      @brief gets the branch from a node pointer (including branch length)
   */ 
  Branch getBranch(nodeptr p) const { return Branch(p->number, p->back->number, p->z[0]); }
  /** 
      @brief gets a nodepointer with specified id 
   */ 
  nodeptr getNode(nat elem) const ; 
  /** 
      @brief extract all branches from the tree (including branch lengths)
   */ 
  std::vector<Branch> extractBranches() const ; 
  /**
     @brief determines the (internal representation of the) tree length
   */ 
  double getTreeLengthExpensive() const;
  

  ///////////////
  // MODIFIERS //
  ///////////////
#if HAVE_PLL != 0 
  void setPartitionList(partitionList *pl) { partitions = *pl; }
  partitionList* getPartitionsPtr()  { return &partitions; }   
#endif  
  /**
     @brief sets the frequencies. Format is important, frequencies must add up to 1.0 
  */ 
  void setFrequencies(const std::vector<double> &values, int model);
  /** 
      @brief sets the parameters. Format is important, last rate must be 1.0  
  */ 
  void setRevMat(const std::vector<double> &values, int model);
  /** 
      @brief sets the alpha for partition "model"
   */ 
  void setAlpha(double alpha,  int model);   
  /** 
      @brief sets a branch. Topology is NOT modified! 
   */ 
  void setBranch(const Branch& b);   
  /** 
      @brief hooks up two nodes with given branch length
   */ 
  void clipNode(nodeptr p, nodeptr q, double z);   
  /** 
      @brief hooks up two nodes with default branch length
   */ 
  void clipNodeDefault(nodeptr p, nodeptr q); 
  /**
     @brief activite parsimony 
   */ 
  void enableParsimony();
  /**
     @brief resets/destroys the topology in the tree 
   */ 
  void unlinkTree();
  /** 
      @brief initializes the tree from a binary file 
  */ 
  void initializeFromByteFile(std::string  byteFileName); 
  /**
     @brief gets a node with given id that is not connected to the tree right now 
   */ 
  nodeptr getUnhookedNode(int number);

  ///////////////
  // observers //
  ///////////////
  std::vector<bool> getExecModel() const ; 
  std::vector<double> getPartitionLnls() const; 
  void setPartitionLnls(const std::vector<double> partitionLnls) ; 
  void setExecModel(const std::vector<bool>  &modelInfo); 

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
  void initRevMat(int model); 	// these functions are not needed any more: directly use the respective setter function     
  void discretizeGamma(int model); 



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
