/**
   @file eval.h
   @brief Functions for likelihood evaluations.  

   New convention: all these functions are responsible for updating
   the chain->likelihood + chain->partitionLnl accordingly. 
*/ 


#ifndef _EVAL_H
#define _EVAL_H

#include "axml.h"
#include "TreeAln.hpp"


class Chain; 



void orientationPointAway(tree *tr, nodeptr p); 
void expensiveVerify(TreeAln& traln); 
void exa_newViewParsimony(TreeAln &traln, nodeptr p); 
void newViewGenericWrapper(TreeAln &traln, nodeptr p, boolean masked); 
void evaluatePartialNoBackup(TreeAln& traln, nodeptr p); 
void evaluateFullNoBackup(TreeAln& traln); 
void evaluateGenericWrapper(TreeAln &traln, nodeptr start, boolean fullTraversal); 
void evaluateOnePartition(TreeAln& traln, nodeptr start, boolean fullTraversal, int model); 

void exa_evaluateParsimony(TreeAln &traln, nodeptr p, boolean fullTraversal, vector<nat> &partitionParsimony); 
void exa_newViewParsimony(TreeAln &traln, nodeptr p); 
#endif
