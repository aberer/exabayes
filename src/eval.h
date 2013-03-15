#ifndef _EVAL_H
#define _EVAL_H

#include "common.h"

void evaluateGenericWrapper(state *chain, nodeptr start, boolean fullTraversal); 
void evaluatePartitions(tree *tr, nodeptr start, boolean fullTraversal, boolean *models); 
void evaluateOnePartition(state *chain, nodeptr start, boolean fullTraversal, int model); 


#endif
