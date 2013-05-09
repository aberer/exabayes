/**
   @file output.h
   
   @brief Functions printing to console or files. 
   
   Also some debug functionality goes in here. 
*/ 

#ifndef  _OUTPUT_H
#define  _OUTPUT_H

#include "TreeAln.hpp"
#include "proposalFunction.h"

class Chain; 
void debug_printTree(TreeAln &traln); 

/* for debugging */
char *Tree2stringNexus(char *treestr, tree *tr , nodeptr p, int perGene ); 

void chainInfoOutput(Chain *chain ); 
void printSample(Chain *chain);
void initializeOutputFiles(Chain *chain);
void finalizeOutputFiles(Chain *chain); 

void chainInfo(Chain *chain); 
void printInfo(Chain *chain, const char *format, ...); 
void debug_printAccRejc(Chain *chain, proposalFunction *pf, boolean accepted) ; 

void debug_printNodeEnvironment(Chain *chain, int nodeID ); 
void debug_checkTreeConsistency(tree *tr); 
void printOrientation(tree *tr, nodeptr p); 


void makeFileNames();

#endif
