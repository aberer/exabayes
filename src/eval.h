#ifndef _EVAL_H
#define _EVAL_H

#include "common.h"
/* #include "bayes.h" */ 


typedef struct 
{
  double likelihood; 		///overall likelihood
  double *partitionLnl; /// likelihood per partition => must be kept up to date 

  double ***vectorsPerPartition; /// stores the lnl arrays  
  char *orientation; /// orientation of the x-vectors for a given partition

} lnlContainer;  /// contains partition and overall
		 /// likelihood. functions in this file have to make
		 /// sure these are always correct. Also contains info to restore lnl computations later 

/* TODO should also contain the posterior */






void evaluateGenericWrapper(struct _state *chain, nodeptr start, boolean fullTraversal); 
void evaluatePartitions(struct _state *chain, nodeptr start, boolean fullTraversal, boolean *models); 
void evaluateOnePartition(struct _state  *chain, nodeptr start, boolean fullTraversal, int model); 


#endif
