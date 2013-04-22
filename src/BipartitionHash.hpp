#ifndef _BIPARTITIONHASH_H

#include "axml.h"
#include "bayes.h"



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


  hashtable* h; 
  nat **bitvectors; 	
  nat numSlots; 
  nat vectorLength; 
  nat numTax; 
  vector<nat> randomHash; 
  vector<string> taxonNames; 
  
}; 


#endif
