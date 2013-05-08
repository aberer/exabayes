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

#include "TreeAln.hpp"


class Chain; 


void evaluateGenericWrapper(Chain *chain, nodeptr start, boolean fullTraversal); 
void evaluateOnePartition(Chain  *chain, nodeptr start, boolean fullTraversal, int model); 
void printAlnTrState(Chain *chain); 
void orientationPointAway(tree *tr, nodeptr p); 
void newViewGenericWrapper(Chain  *chain, nodeptr p, boolean masked); 

void evaluateFullNoBackup(Chain *chain); 
void evaluatePartialNoBackup(Chain *chain, nodeptr p); 

void expensiveVerify(Chain *chain); 
nat exa_evaluateParsimony(TreeAln &traln, nodeptr p, boolean fullTraversal ); 
void exa_newViewParsimony(TreeAln &traln, nodeptr p); 
#endif
