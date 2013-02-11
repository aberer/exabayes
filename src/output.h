#ifndef  _OUTPUT_H
#define  _OUTPUT_H

void makeFileNames(void); 
void chainInfoOutput(state *curstate, int sum_radius_accept, int sum_radius_reject ); 
void printSample(state *curstate);
void initializeOutputFiles(); 
void finalizeOutputFiles(state *curstate); 

#endif







