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

#if HAVE_PLL == 0  
  void initializePartitionsExaml(TreeAln &traln, std::ifstream &byteStream); 
  void initializeTreeExaML(TreeAln &traln, std::string byteFileName ); 
#else 
  void initializeTreePLL(TreeAln &traln, std::string byteFileName);
  void initializePartitionsPLL(TreeAln &traln, std::string byteFileName);
#endif

  void readPartitions(TreeAln & traln, std::ifstream &byteStream); 

  void initializeParsimonyVectors(TreeAln &traln, std::ifstream& byteStream); 
  void unifiedInitializePartitions(TreeAln &traln, std::string byteFileName); 
  void parseMagicNumber(std::ifstream& in); 

private:
  void unifiedModelInit(TreeAln &traln); 

#if HAVE_PLL == 0
  analdef adef; 

#else 
#endif
  
  unsigned int mask32[32]; 
}; 



#endif
