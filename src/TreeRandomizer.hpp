
#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H
#include <memory>

#include "axml.h" 
#include "Randomness.hpp"
#include "TreeAln.hpp"

class TreeRandomizer
{
public: 
  TreeRandomizer(int seed);
  void randomizeTree( TreeAlnPtr &traln);
  
private: 
  Randomness rand; 

  int markBranches(TreeAlnPtr &traln, nodeptr *branches, nodeptr p, int *counter, int numsp); 
  void insertTaxon ( TreeAlnPtr &traln, nodeptr p, nodeptr q); 
  void buildSimpleTreeRandom (TreeAlnPtr &traln, int ip, int iq, int ir); 
  nodeptr buildNewTip (TreeAlnPtr &traln, nodeptr p); 
  


}; 


#endif
