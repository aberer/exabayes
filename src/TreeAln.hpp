/** 
    @brief mostly wraps the legacy tr and partition list     
*/


#ifndef _TREEALN_H
#define _TREEALN_H

#include <vector>
#include <iostream>
#include <array>
#include <memory>

#include "axml.h"
#include "GlobalVariables.hpp"

#include "BranchFwd.hpp"
#include "FlagType.hpp"
#include "RunModes.hpp"

class AbstractPrior; 
class Randomness; 
class TreePrinter; 
class AbstractParameter; 

std::ostream& operator<<(std::ostream& out, pInfo& rhs); 
nat numStateToNumInTriangleMatrix(int numStates) ; 


#if HAVE_PLL == 0 
struct partitionList; 
#endif

class TreeAln
{
public: 
   /////////////////
   // life cycle  //
   /////////////////
  explicit TreeAln();
  ~TreeAln();
  /** 
      @brief copies the entire model from the rhs.
      Mostly for debugging purposes. 
   */ 
  void copyModel(const TreeAln& rhs)  ; 
  /**
     @brief copies the entire state from the rhs to this tree/alignment.
     Important: likelihood arrays are not concerned
  */ 
  TreeAln& operator=( TreeAln &rhs) = delete ;     
  TreeAln(const TreeAln &tmp) = delete ; 
  
  /** 
      @brief indicates whether this TreeAln is equal to rhs (with
      regard to all model parameters)
   */   
  bool operator==(const TreeAln& rhs); 


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
      @brief frees all likelihood arrays 
   */   
  void clearMemory(); 
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
  bool isTipNode(nat node) const { return isTip(node, getTr()->mxtips) ; }
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
  BranchLength getBranch(const BranchPlain &branch,  const AbstractParameter *param) const; 
  BranchLengths getBranch(const BranchPlain& branch, const std::vector<AbstractParameter*> &params) const ; 

  BranchLength getBranch(const nodeptr &branch,  const AbstractParameter *param) const; 
  BranchLengths getBranch(const nodeptr &branch, const std::vector<AbstractParameter*> &params) const ; 
  /** 
      @brief gets a nodepointer with specified id 
   */ 
  nodeptr getNode(nat elem) const ; 
  /** 
      @brief extract all branches from the tree (including branch lengths)
   */ 
  // template<typename RESULT>
  std::vector<BranchLength> extractBranches( const AbstractParameter* param ) const; 
  std::vector<BranchLengths> extractBranches( const std::vector<AbstractParameter*> &param ) const;  
  std::vector<BranchPlain> extractBranches() const; 
  /** 
      @brief gets the number of inner nodes in the tree 
   */ 
  nat getNumberOfInnerNodes() const { return getNumberOfNodes()  - getNumberOfTaxa()  ;   } 
  
  /** 
      @brief gets the mean substitution rate overall specified partitions
   */ 
  double getMeanSubstitutionRate(const std::vector<nat> &partitions) const ;
  /** 
      @brief unlinks a node 
   */ 
  void unlinkNode(nodeptr p); 
  /** 
      @brief gets the three nodes adjacent to the given node  
   */ 
  std::vector<nat> getNeighborsOfNode( nat node ) const ; 
  /** 
      @brief prunes the node from the tree 
   */ 
  void detachNode(nodeptr p); 
  ///////////////
  // MODIFIERS //
  ///////////////
#if HAVE_PLL != 0 
  void setPartitionList(partitionList *pl) { partitions = *pl; }
  partitionList* getPartitionsPtr()  { return &partitions; }   
#endif  
  BranchPlain getAnyBranch() const ;//  {return BranchPlain(tr.nodep[1]->number, tr.nodep[1]->back->number); } 
  
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
  void setBranch(const BranchLengths& b, const std::vector<AbstractParameter*> params);   
  void setBranch(const BranchLength& branch, const AbstractParameter* param); 
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
      @brief gets the maximum length of paths below this branch 
   */ 
  nat getDepth(const BranchPlain  &b) const ; 
  /** 
      @brief gets the longest path 
   */ 
  std::vector<nat> getLongestPathBelowBranch(const BranchPlain &b) const ; 
  /** 
      @brief initializes the tree from a binary file.  It guarantees
      that the tree is in a usable state afterwards (and thus, this
      method may me rather expensive)
  */ 
  void initializeFromByteFile(std::string byteFileName, RunModes mode); 
  void initializeNumBranchRelated(); 
  // nat readNUM_BRANCHES(std::string byteFileName); 
  
  /**
     @brief gets a node with given id that is not connected to the tree right now 
   */ 
  nodeptr getUnhookedNode(int number);
  ///////////////
  // observers //
  ///////////////
  /** 
      @brief gets the branches that are below one branch 

      This function is handy for traversing the tree and relying less
      on raw pointers. For traversing, it is often necessary to invert
      the resulting branch.      
   */
  std::pair<BranchPlain,BranchPlain> getDescendents(const BranchPlain &b) const; 

  std::vector<bool> getExecModel() const ; 
  std::vector<double> getPartitionLnls() const; 
  void setPartitionLnls(const std::vector<double> partitionLnls) ; 
  void setExecModel(const std::vector<bool>  &modelInfo); 
  /** 
      @brief gets a list of of sequence names  
   */ 
  std::vector<std::string> getNameMap() const; 

  static const double zZero;   
  static const double initBL;  	// init values 
  static const double problematicBL; 

  nat getNumberOfAssignedSites(nat model) const ; 

  bool revMatIsImmutable(int model) const; 

  RunModes getMode() const { return mode; }

private: 			// METHODS  
#if HAVE_PLL != 0

#endif  

  // template<class BRANCH> void extractHelper(nodeptr p, std::vector<BRANCH> &result, bool isStart) const; 

  void initRevMat(int model); 	// these functions are not needed any more: directly use the respective setter function     
  void discretizeGamma(int model); 

private: 			// ATTRIBUTES 
#if HAVE_PLL != 0 
  partitionList partitions; 
#endif
  tree tr;		// TODO replace with an object for cleanup   
  bool parsimonyEnabled;   
  RunModes mode; 
  
  friend std::ostream& operator<< (std::ostream& out,  const TreeAln&  traln);
};  


#endif
