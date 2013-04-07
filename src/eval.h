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

typedef struct 
{
  double ***vectorsPerPartition; /// stores the lnl arrays  
  int *orientation;  /// indicates the node towards a specific node is directed  
  
  nat **partitionScaler;  /// for each partition, for each node

} lnlContainer;  /// contains partition and overall
		 /// likelihood. functions in this file have to make
		 /// sure these are always correct. Also contains info to restore lnl computations later 

/* TODO should also contain the posterior */

struct _state; 


void evaluateGenericWrapper(struct _state *chain, nodeptr start, boolean fullTraversal); 
void evaluatePartitions(struct _state *chain, nodeptr start, boolean fullTraversal, boolean *models); 
void evaluateOnePartition(struct _state  *chain, nodeptr start, boolean fullTraversal, int model); 

void saveOrientation(struct _state *chain); 
void loadOrientation(struct _state *chain); 
void saveArray(struct _state *chain, int model); 
void loadArray(struct _state *chain, int model); 
void restoreAlignAndTreeState(struct _state *chain); 
void printAlnTrState(struct _state *chain); 
void saveAlignAndTreeState(struct _state *chain); 
void orientationPointAway(tree *tr, nodeptr p); 
#endif
