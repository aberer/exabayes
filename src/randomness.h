#ifndef _RANDOMNESS_H
#define _RANDOMNESS_H

#include "rng.h"


/* todo replace that once */
typedef struct _rngState
{
  randKey_t key; 
  randCtr_t ctr; 
} rngState; 


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

branch drawBranchUniform(state *chain); 
#endif
