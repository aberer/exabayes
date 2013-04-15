/**
   @file convergence.h
   
   @brief Convergence diagnostics.
*/ 


#ifndef _CONVERGENCE_H 
#define _CONVERGENCE_H 

void addBipartitionsToHash(tree *tr, state *chain);
boolean convergenceDiagnostic(state *allChains); 
void initializeConvergenceStructs(tree *tr, state *chain); 

#endif
