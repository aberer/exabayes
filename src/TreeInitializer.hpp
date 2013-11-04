/**
   @file TreeInitializer.hpp
   
   @brief basic functionality to set up a tree. 

   The initialization of raxml's tree-struct is somewhat involved.
 */ 


#ifndef _TREE_INITIALIZER_HPP
#define _TREE_INITIALIZER_HPP

#include <memory>
#include <string>
#include "InitializationResource.hpp"
#include "common.h"
#include "RunModes.hpp"
#include "axml.h"

class TreeAln; 

class TreeInitializer
{
public: 
  TreeInitializer(std::unique_ptr<InitializationResource> initRes); 

  void initializeWithAlignmentInfo(TreeAln &traln, RunModes flags); 
  static void setupTheTree(tree &tr); 
  static void initializeBranchLengths(tree &tr, nat numPart, nat numTax ); 

  void determinePartitionContribution(TreeAln &traln); 

private:			// METHODS
  void unifiedInitializePartitions(TreeAln &traln); 
  /** 
      @brief sets the weights and yVectors 
   */ 
  void distributeData(TreeAln &traln);
#if HAVE_PLL == 0  
  void initializePartitionsExaml(TreeAln &traln); 
  void initializeTreeExaML(TreeAln &traln ); 
#else 
  void initializeTreePLL(TreeAln &traln);
  void initializePartitionsPLL(TreeAln &traln);
#endif
  void initializeParsimonyVectors(TreeAln &traln); 
  void copyWeightsAndAln(TreeAln &traln );
  void parseMagicNumber(std::ifstream& in); 
  void unifiedModelInit(TreeAln &traln); 

private:  			// ATTRIBUTES
#if HAVE_PLL == 0
  analdef adef; 
#endif
  unsigned int mask32[32]; 

  std::unique_ptr<InitializationResource> _initResPtr; 
}; 

#endif
