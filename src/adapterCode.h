#ifndef _ADAPTER_CODE_H
#define _ADAPTER_CODE_H

#include "config.h"


/* shared  */
void exa_newViewGeneric(tree *tr, nodeptr p, boolean masked); 
void exa_initReversibleGTR( tree *tr,int model); 
void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 
void exa_evaluateGeneric(tree *tr, nodeptr start, boolean fullTraversal); 


int getNumberOfPartitions(tree *tr) ; 
pInfo* getPartition(tree *tr, int num); 
boolean hasPergeneBL(tree *tr); 

double getFracChange(tree *tr, int num); 
void setFracChange(tree *tr, int num, double value); 

double getPLH(tree *tr, int num) ; 
void setPLH(tree *tr, int num, double value); 


void setExecModel(tree *tr, int num,boolean value); 
boolean getExecModel(tree *tr, int num); 


#endif
