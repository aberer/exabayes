#ifndef  _OUTPUT_H
#define  _OUTPUT_H

/* for debugging */
char *Tree2stringNexus(char *treestr, tree *tr, nodeptr p, int perGene ); 

void makeFileNames(void); 
void chainInfoOutput(state *curstate ); 
void printSample(state *curstate);
void initializeOutputFiles(); 
void finalizeOutputFiles(state *curstate); 

#endif







