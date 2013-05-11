/**
   @file adapters.h

   @brief All code that ensures that ExaBayes can be compiled with
   either the PLL or the ExaML code base

*/ 


#ifndef _ADAPTER_CODE_H
#define _ADAPTER_CODE_H

void exa_newViewGeneric(TreeAln& traln, nodeptr p, boolean masked); 
void exa_hookupDefault(tree *tr, nodeptr p, nodeptr q); 
void exa_evaluateGeneric(TreeAln &traln, nodeptr start, boolean fullTraversal); 
bool isOutputProcess(); 
#endif
