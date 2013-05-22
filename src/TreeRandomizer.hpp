
#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H
#include <memory>

#include "axml.h" 
#include "Randomness.hpp"


class TreeRandomizer
{
public: 
  TreeRandomizer(int seed, shared_ptr<TreeAln> traln);
  shared_ptr<TreeAln> getTr(){return traln; }
  void randomizeTree();
  
private: 
  Randomness rand; 
  shared_ptr<TreeAln> traln; 

  int markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp); 
  void insertTaxon ( nodeptr p, nodeptr q); 
  void buildSimpleTreeRandom (int ip, int iq, int ir); 
  nodeptr buildNewTip (nodeptr p); 
  


}; 


#endif
