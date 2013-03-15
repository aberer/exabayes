#ifndef _ADAPTER_CODE_H
#define _ADAPTER_CODE_H

#include "config.h"


/* shared  */
void exa_newViewGeneric(state *chain, nodeptr p, boolean masked); 
void exa_initReversibleGTR( tree *tr,int model); 
void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 
void exa_evaluateGeneric(tree *tr, nodeptr start, boolean fullTraversal); 

int getNumBranches(tree *tr); 
void setNumbranches(tree *tr, int num); 

void setNumberOfPartitions(tree *tr, int num); 
pInfo* getPartition(state *chain, int num); 
int getNumberOfPartitions(tree *tr) ; 

boolean hasPergeneBL(tree *tr); 

void setPLH(state *chain, int num, double value); 
double getPLH(state *chain, int num); 

void setExecModel(state *chain, int num,boolean value); 
boolean getExecModel(state *chain, int num); 

double getPcontr(state *chain, int num); 
double getFracChange(state *chain, int num); 
void setFracChange(state *chain, int num, double value); 



#endif
