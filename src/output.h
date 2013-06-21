/**
   @file output.h
   
   @brief Functions printing to console or files. 
   
   Also some debug functionality goes in here. 
*/ 

#ifndef  _OUTPUT_H
#define  _OUTPUT_H

#include "TreeAln.hpp"

class Chain; 
void debug_printTree(TreeAln &traln); 

void chainInfoOutput(Chain *chain ); 
void printSample(Chain *chain);
void initializeOutputFiles(Chain *chain);
void finalizeOutputFiles(Chain *chain); 

void chainInfo(Chain *chain); 
void printInfo(Chain *chain, const char *format, ...); 

void debug_checkTreeConsistency(const TreeAln& traln); 
void printOrientation(tree *tr, nodeptr p); 

#endif
