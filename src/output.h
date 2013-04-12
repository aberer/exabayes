/**
   @file output.h
   
   @brief Functions printing to console or files. 
   
   Also some debug functionality goes in here. 
*/ 

#ifndef  _OUTPUT_H
#define  _OUTPUT_H

#ifdef __cplusplus
extern "C"{
#endif


void debug_printTree(state *chain); 

/* for debugging */
char *Tree2stringNexus(char *treestr, state *chain , nodeptr p, int perGene ); 

/* void makeFileNames();  */
void chainInfoOutput(state *chain ); 
void printSample(state *chain);
void initializeOutputFiles(state *chain);
void finalizeOutputFiles(state *chain); 

void chainInfo(state *chain); 
void printInfo(state *chain, const char *format, ...); 
void debug_printAccRejc(state *chain, proposalFunction *pf, boolean accepted) ; 

void debug_printNodeEnvironment(state *chain, int nodeID ); 
void debug_checkTreeConsistency(state *chain); 

void printOrientation(tree *tr, nodeptr p); 

#ifdef __cplusplus
}
#endif
#endif
