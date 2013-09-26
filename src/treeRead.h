/**
   @file treeRead.h 
   @brief functions for reading trees 
   
   Mostly taken over from the axml-variants.  
 
*/ 

#ifndef _TREE_READ_H
#define _TREE_READ_H

/* boolean readTreeWithOrWithoutBL(tree *tr, std::string treeString);  */
/* void myTreeReadLen(FILE *fp, tree *tr, boolean hasBL); */
boolean readTreeWithOrWithoutBL(tree *tr, std::string treeString); 
void myTreeReadLen(std::string treeString , tree *tr, boolean hasBL); 

#endif
