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

typedef struct _state state; 

void evaluateGenericWrapper(struct _state *chain, nodeptr start, boolean fullTraversal); 
void evaluateOnePartition(struct _state  *chain, nodeptr start, boolean fullTraversal, int model); 
void printAlnTrState(struct _state *chain); 
void orientationPointAway(tree *tr, nodeptr p); 
void newViewGenericWrapper(struct _state  *chain, nodeptr p, boolean masked); 

void evaluateFullNoBackup(state *chain); 
void evaluatePartialNoBackup(state *chain, nodeptr p); 

void expensiveVerify(state *chain); 
#endif
