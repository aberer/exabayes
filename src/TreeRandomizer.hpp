
#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H
#include <memory>

#include "axml.h" 
#include "Randomness.hpp"
#include "TreeAln.hpp"

class TreeRandomizer
{
public: 
  TreeRandomizer(randCtr_t seed);
  void randomizeTree(TreeAln &traln);
  
private: 
  Randomness rand; 
}; 


#endif
