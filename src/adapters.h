/**
   @file adapters.h

   @brief All code that ensures that ExaBayes can be compiled with
   either the PLL or the ExaML code base

*/ 


#ifndef _ADAPTER_CODE_H
#define _ADAPTER_CODE_H

#include "config.h"


/* shared  */
void exa_evaluateGeneric(Chain *chain, nodeptr start, boolean fullTraversal); 
void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 
void exa_newViewGeneric(Chain *chain, nodeptr p, boolean masked); 
bool isOutputProcess(); 
#endif
