#ifndef _TREE_RANDOMIZER_H
#define _TREE_RANDOMIZER_H

#include "axml.h" 
#include "Randomness.hpp"
#include "TreeAln.hpp"

class TreeRandomizer
{
public: 

  /**
     @brief creates a random tree 
   */   
  static void randomizeTree( TreeAln &traln, Randomness &rand);


  static Branch drawBranchUniform_helper(const TreeAln &traln, Randomness &rand , nat curNumTax)  ; 
  static Branch drawBranchUniform(const TreeAln & traln, Randomness &rand )  ; 
  static Branch drawInnerBranchUniform( const TreeAln& traln, Randomness &rand)  ; 
  static Branch drawBranchWithInnerNode(const TreeAln& traln, Randomness &rand)  ; 
  static nat drawInnerNode(const TreeAln& traln, Randomness &rand )  ; 

  /**
     @brief creates a parsimony tree and applies it to the specfilied tree.      
     the routine does not enforce treatment of branch lengths.  
   */ 
  static void createParsimonyTree(TreeAln &traln, Randomness& rand); 

}; 


#endif
