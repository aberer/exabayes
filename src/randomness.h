#ifndef _RANDOMNESS_H
#define _RANDOMNESS_H

void initLocalRng(state *theChain); 

int drawGlobalRandIntBound(int upperBound); 
randCtr_t drawGlobalRandInt(); 
double drawGlobalDouble01();

int drawRandInt(state *chain, int upperBound); 
double drawRandDouble01(state *chain);
double drawRandExp(state *chain, double lambda);
double drawRandBiUnif(state *chain, double x);
int drawSampleProportionally(state *chain,  double *weights, int numWeight ); 
void drawPermutation(state *chain, int* perm, int n); 
#endif








