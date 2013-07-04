#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H

#include "axml.h" 
#include "Randomness.hpp"
#include "TreeAln.hpp"

class TreeRandomizer
{
public: 
  TreeRandomizer(randCtr_t seed);

  /**
     @brief creates a random tree 
   */   
  void randomizeTree(TreeAln &traln);


  /**
     @brief creates a parsimony tree and applies it to the specfilied tree.      
     the routine does not enforce treatment of branch lengths.  
   */ 
  void createParsimonyTree(TreeAln &traln); 
  
private: 
  Randomness rand; 
}; 


#endif
