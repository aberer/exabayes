/**
   @file eval.h
   @brief Functions for likelihood evaluations.  

   New convention: all these functions are responsible for updating
   the chain->likelihood + chain->partitionLnl accordingly. 
*/ 


#ifndef _EVAL_H
#define _EVAL_H

#include "axml.h"
#include "common.h"

typedef struct  _savedArray
{
  int node; 			/* node this array was associated to */
  struct _savedArray *next; 	
} savedArray; 


/* typedef struct  */
/* { */
  /* double ***vectorsPerPartition; /// stores the lnl arrays   */
  /* int *orientation;  /// indicates the node towards a specific node is directed   */
  
  /* nat **partitionScaler;  /// for each partition, for each node */
  
  /* savedArray *savedArrayList; 	/\* contains information, if a switch was done with a respective node  *\/ */
  /* int partitionEvaluated; 	/\* which parition was evaluated, if there has been only one  *\/ */
  /* int numberArrays; 		/\* number of arrays taht have been exchanged (i.e., for how many parttions)   *\/ */
  
/* } lnlContainer;  /// contains partition and overall */
		 /// likelihood. functions in this file have to make
		 /// sure these are always correct. Also contains info to restore lnl computations later 

/* TODO should also contain the posterior */


typedef struct _state state; 


void evaluateGenericWrapper(struct _state *chain, nodeptr start, boolean fullTraversal); 
void evaluateOnePartition(struct _state  *chain, nodeptr start, boolean fullTraversal, int model); 
void printAlnTrState(struct _state *chain); 
void orientationPointAway(tree *tr, nodeptr p); 
void newViewGenericWrapper(struct _state  *chain, nodeptr p, boolean masked); 
void evaluateFullNoBackup(state *chain); 

#endif
