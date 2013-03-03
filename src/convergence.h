#ifndef _CONVERGENCE_H 
#define _CONVERGENCE_H 

void addBipartitionsToHash(tree *tr, state *chain);
boolean convergenceDiagnostic(state *allChains, int numChains); 
void initializeConvergenceStructs(tree *tr, state *chain); 
#endif
