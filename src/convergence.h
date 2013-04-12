/**
   @file convergence.h
   
   @brief Convergence diagnostics.
*/ 


#ifndef _CONVERGENCE_H 
#define _CONVERGENCE_H 



/* #ifdef __cplusplus */
/* extern "C"{ */
/* #endif */

void addBipartitionsToHash(tree *tr, state *chain);
boolean convergenceDiagnostic(state *allChains); 
void initializeConvergenceStructs(tree *tr, state *chain); 


/* #ifdef __cplusplus */
/* } */
/* #endif */
#endif
