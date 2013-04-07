/**
   @file adapters.h

   @brief All code that ensures that ExaBayes can be compiled with
   either the PLL or the ExaML code base

*/ 


#ifndef _ADAPTER_CODE_H
#define _ADAPTER_CODE_H

#include "config.h"


#if HAVE_PLL == 1 
void initializeTree(tree *tr, partitionList *partitions, analdef *adef); 
#else  
void initializeTree(tree *tr, analdef *adef); 
#endif


/* shared  */
void exa_evaluateGeneric(state *chain, nodeptr start, boolean fullTraversal); 
void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 
void exa_newViewGeneric(state *chain, nodeptr p, boolean masked); 
void exa_initReversibleGTR( state *chain,int model); 

int getNumBranches(tree *tr);
pInfo* getPartition(state *chain, int num); 
int getNumberOfPartitions(tree *tr) ; 
boolean hasPergeneBL(tree *tr); 
void setPLH(state *chain, int num, double value); 
double getPLH(state *chain, int num); 
void setExecModel(state *chain, int num,boolean value); 
boolean getExecModel(state *chain, int num); 
double getPcontr(state *chain, int num); 
/* double getFracChange(state *chain, int num);  */

double** getXPtr(state *chain, int model); 
#endif
