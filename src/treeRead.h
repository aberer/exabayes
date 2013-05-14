/**
   @file treeRead.h 
   @brief functions for reading trees 
   
   Mostly taken over from the axml-variants.  
 
*/ 

#ifndef _TREE_READ_H
#define _TREE_READ_H

void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 
boolean readTreeWithOrWithoutBL(tree *tr, FILE *treeFile); 
void myTreeReadLen(FILE *fp, tree *tr, boolean hasBL); 
void traverseInitCorrect(nodeptr p, int *count, TreeAln *traln ); 

#endif
