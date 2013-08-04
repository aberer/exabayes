#ifndef _TREEALN_H
#define _TREEALN_H

#include <vector>
#include <iostream>
#include <memory>

#include "priors/AbstractPrior.hpp"
#include "Branch.hpp"
#include "axml.h"
#include "Randomness.hpp"
#include "GlobalVariables.hpp"
#include "TreePrinter.hpp"

std::ostream& operator<<(std::ostream& out, pInfo& rhs); 

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
  TreeAln& operator=( TreeAln &rhs) = delete ; 
  TreeAln(const TreeAln &tmp) = delete ; 


  /////////////////////////////////////
  //           OBSERVERS             //
  /////////////////////////////////////
  /**
     @brief get the raxml representation of the partition  
   */ 
  pInfo* getPartition(nat model) const;
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
      @notice do NOT confuse with getNumberOfBranches 
      @todo remodel this function / deal with the general problem 
   */ 
  nat getNumBranches() const; 
  /** 
      @brief get the number of branches in the tree (not counting per-partition branch lengths)
   */ 
  nat getNumberOfBranches() const {return getNumberOfNodes() - 1 ; }
  /** 
      @brief get the number of partitions in the alignment
   */ 
  nat getNumberOfPartitions() const;   
  /** 
      @brief get the number of taxa in the tree 
   */ 
  nat getNumberOfTaxa() const {return getTr()->mxtips; }
  /** 
      @brief get the number of nodes in the unrooted tree 
   */ 
  nat getNumberOfNodes() const { nat numTax = getNumberOfTaxa(); return 2 * numTax - 2 ;  } // excluding the virtual root 
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
      
      @notice if only one parameter is provided, the resulting branch
      will ONLY contain the branch length for this one parameter
      (bogus parameters for other branch lengths, if you have
      per-partition branch lengths). For setting this branch, the same
      parameter must be used (and only this one).
      
      If you call this function with an existing Branch b, this
      function is usefull to get the actual branch length.
   */ 
  Branch getBranch(nodeptr p, const std::vector<AbstractParameter*> &param) const; 
  Branch getBranch(const Branch& b, AbstractParameter* param); 
  Branch getBranch(nodeptr p, const AbstractParameter* param) const; 
  Branch getBranch(const Branch &b , const std::vector<AbstractParameter*> params) const ; 
  /** 
      @brief gets a nodepointer with specified id 
   */ 
  nodeptr getNode(nat elem) const ; 
  /** 
      @brief extract all branches from the tree (including branch lengths)
   */ 
  std::vector<Branch> extractBranches(const std::vector<AbstractParameter*> &params) const ; 
  std::vector<Branch> extractBranches( AbstractParameter* params) const ; 
  /**
     @brief determines the (internal representation of the) tree length
   */ 
  std::vector<double> getTreeLengthsExpensive(const std::vector<AbstractParameter*> &params) const;  
  /** 
      @brief gets the number of inner nodes in the tree 
   */ 
  nat getNumberOfInnerNodes() const { return getNumberOfNodes()  - getNumberOfTaxa()  ;   } 
  
  /** 
      @brief gets the mean substitution rate overall specified partitions
   */ 
  double getMeanSubstitutionRate(const std::vector<nat> &partitions) const ;

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
  void setBranch(const Branch& b, const std::vector<AbstractParameter*> params);   
  void setBranch(const Branch& branch, AbstractParameter* param); 
  /** 
      @brief hooks up two nodes. Branch lengths must be set
      separately.
   */ 
  void clipNode(nodeptr p, nodeptr q);   
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
      @brief initializes the tree from a binary file.  It guarantees
      that the tree is in a usable state afterwards (and thus, this
      method may me rather expensive)
  */ 
  void initializeFromByteFile(std::string  byteFileName); 
  /**
     @brief gets a node with given id that is not connected to the tree right now 
   */ 
  nodeptr getUnhookedNode(int number);
  /** 
      @brief copies the entire model from the rhs.
      Mostly for debugging purposes. 
   */ 
  void copyModel(const TreeAln& rhs); 
  ///////////////
  // observers //
  ///////////////
  /** 
      @brief gets the branches that are below one branch 
   */
  std::pair<Branch,Branch> getDescendingBranches(const Branch &b) const; 
  std::vector<bool> getExecModel() const ; 
  std::vector<double> getPartitionLnls() const; 
  void setPartitionLnls(const std::vector<double> partitionLnls) ; 
  void setExecModel(const std::vector<bool>  &modelInfo); 

  static const double zZero;   
  static const double initBL;  	// init values 
  static const double problematicBL; 

  // DEBUG 
  void printArrayStart(std::ostream &out, nat length = 2 ); 

private: 			// METHODS  
#if HAVE_PLL != 0
  void initializeTreePLL(std::string byteFileName);
  void initializePartitionsPLL(std::string byteFileName, double ***empFreq, bool multiBranch);
#endif  

  double getTreeLengthHelper(nodeptr p) const;
  void extractHelper( nodeptr p , std::vector<Branch> &result, bool isStart, const std::vector<AbstractParameter*> &params) const ; 
  void extractHelper( nodeptr p , std::vector<Branch> &result, bool isStart, AbstractParameter* param ) const ; 
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
  bool revMatIsImmutable(int model) const ; 

  bool usingPerSiteRates()  { return tr.rateHetModel == GAMMA;  } 
  void enablePerSiteRates() {  tr.rateHetModel = CAT; } 

};  


#endif
