#ifndef _TREE_INITIALIZER_HPP
#define _TREE_INITIALIZER_HPP

#include <string>
#include "common.h"
#include "axml.h"

class TreeAln; 

class TreeInitializer
{
public: 
  TreeInitializer(); 

  void unifiedInitializePartitions(TreeAln &traln, std::string byteFileName); 

private:			// METHODS
#if HAVE_PLL == 0  
  void initializePartitionsExaml(TreeAln &traln, std::ifstream &byteStream); 
  void initializeTreeExaML(TreeAln &traln, std::string byteFileName ); 
#else 
  void initializeTreePLL(TreeAln &traln, std::string byteFileName);
  void initializePartitionsPLL(TreeAln &traln, std::string byteFileName);
#endif
  void initializeParsimonyVectors(TreeAln &traln, std::ifstream& byteStream); 
  void readPartitions(TreeAln & traln, std::ifstream &byteStream); 
  void parseMagicNumber(std::ifstream& in); 
  void unifiedModelInit(TreeAln &traln); 

private:  			// ATTRIBUTES
#if HAVE_PLL == 0
  analdef adef; 
#endif
  unsigned int mask32[32]; 
}; 

#endif
