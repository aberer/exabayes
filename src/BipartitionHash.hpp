#ifndef _BIPARTITIONHASH_H

#include <vector>

#include "axml.h"
#include "TreeAln.hpp"

class BipartitionHash
{
public: 
  BipartitionHash(int numTax, int numRuns);
  ~BipartitionHash();

  void addBipartitionsToHash(TreeAln &taln, nat chainid); 
  double averageDeviationOfSplitFrequencies(double ignoreFreq); 

  
private: 
  void newviewBipartitions( nodeptr p ); 
  void getxnodeBips (nodeptr p); 
  void insertAndCount(TreeAln &traln, nat *bitVector, hashNumberType position, int chainId); 
  void extractBipartitions(TreeAln &traln, nodeptr p, int *cnt, int chainId); 
  void resetBitVectors();
  void initializeRandomHash(); 
  void printBv(nat *bv); 


  hashtable* h; 
  nat **bitvectors; 	
  nat numSlots; 
  nat vectorLength; 
  nat numTax; 
  std::vector<nat> randomHash; 
  std::vector<std::string> taxonNames; 
  
}; 


#endif
