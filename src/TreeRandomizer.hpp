
#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H

#include "axml.h" 
#include "bayes.h"
#include "Randomness.hpp"


class TreeRandomizer
{
public: 
  TreeRandomizer(int seed, TreeAln* traln);
  TreeAln *getTr(){return traln; }
  void randomizeTree();
  
private: 
  Randomness rand; 
  TreeAln *traln; 

  int markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp); 
  void insertTaxon ( nodeptr p, nodeptr q); 
  void buildSimpleTreeRandom (int ip, int iq, int ir); 
  nodeptr buildNewTip (nodeptr p); 
  


}; 


#endif
