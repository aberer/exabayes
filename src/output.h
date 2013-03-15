#ifndef  _OUTPUT_H
#define  _OUTPUT_H

/* for debugging */
char *Tree2stringNexus(char *treestr, state *chain , nodeptr p, int perGene ); 

void makeFileNames(void); 
void chainInfoOutput(state *curstate ); 
void printSample(state *curstate);
void initializeOutputFiles(); 
void finalizeOutputFiles(state *curstate); 

void chainInfo(state *chain); 

#endif
